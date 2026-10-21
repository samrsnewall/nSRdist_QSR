function [nSRcounts, agediffs] = oneCoreRSR(corename, dataLoc, scenarios, LabIDs, incDepths, excLabIDs, excDepths, scenario_meanSR, ageModes, S, minAgeDiff, sampleOrder)
% This function calculates nSR for one core by random sampling of the
% radiocarbon ages (the RSR methods). For each of S.numruns runs it picks one
% age-depth scenario at random, samples one age per date from its calibrated
% PDF in the order sampleOrder while forcing ages to increase with depth,
% drops any date within minAgeDiff yr of the last kept date, and works out the
% nSR of every remaining pair of dates. minAgeDiff = 0, 500, 1000 give
% RSR0, RSR500, RSR1000. Raw 14C ages and errors are in kyr; the
% calibrated age grid and everything sampled from it are in yr; depths in cm;
% SR in cm/kyr.
%
% scenarios holds one string array of LabIDs per valid scenario ("invalid"
% marks a rejected one); scenario_meanSR (cm/kyr) is used to normalise when
% S.normWithRunAve is false; ageModes gives the mode age of each date per
% scenario for S.useModes. nSRcounts is the 3-row nSR matrix (nSR, dx_cm,
% dt_yr) with a NaN column starting each run; agediffs is every dt_yr from
% every run. sampleOrder is "shallow2deep" (the manuscript method),
% "deep2shallow", "random" or "rejection" (joint rejection sampling, with
% the minAgeDiff restriction applied shallow to deep afterwards).
%
% S fields used: numruns, useModes, normWithRunAve, pdfMinVal, rejectionBatchSize,
% rejectionMaxBatches

%% Check whether the core has previously been rejected
% Stop here if the core has been previously rejected
if isempty(scenarios)
    % Core was rejected previously
    nSRcounts = [];
    agediffs  = [];
    return
end

%% Read in radiocarbon data
% Read in the radiocarbon ages, depths, errors and LabIDs for this core
% (from the World Atlas NetCDF file or the Lin2014 text file)
if dataLoc == "WA"
    [age_kyr, depth_cm, ageErr_kyr, label] = getDataWA(corename, S);
elseif dataLoc == "Lin2014"
    [age_kyr, depth_cm, ageErr_kyr, label] = getDatatxt(corename, S);
end

%% Filter ages
% Filter the ages by material, manual exclusions and calibration curve range
[age_kyr, depth_cm, ageErr_kyr, label, noIncludeInfo, tooFewAges] = filtering(age_kyr, ...
    depth_cm, ageErr_kyr, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

% Stop here if filtering of ages leads to core rejection
if noIncludeInfo == 1 || tooFewAges == 1
    % Core was rejected after filtering through ages
    nSRcounts = [];
    agediffs  = [];
    return
end

%% Create fake labels if needed
% If there are no LabIDs, create temporary LabIDs so that scenarios can be made
if sum(contains(label, "NaN")) ~= 0
    for i = find(contains(label, "NaN"))
        label(i) = "FakeLabel" + num2str(i);
    end
end

% If any label is duplicated, replace every label with a numbered fake one,
% exactly as oneCoreScenarios did when it built the scenarios, so that the
% labels here match the scenario labels one to one. The dates are in the
% same depth order in both, because both come from the same filtering call
[~, uniqueIdx] = unique(label);
if length(uniqueIdx) ~= length(label)
    for iLabel = 1:length(label)
        label(iLabel) = "xFakeLabel" + num2str(iLabel);
    end
end

%% Check which scenarios to analyse
% Initialise every scenario as valid
validScenariosLog = true(length(scenarios), 1);

% Check which scenarios are invalid (due to reversals etc)
for i_sce = 1:length(scenarios)
    if scenarios{i_sce} == "invalid"
        validScenariosLog(i_sce) = false;
    end
end

%% Loop through scenarios using random sampling SR calculations
% Calibrate the ages (one probability column per date, on the calAge_yr grid)
[ageProbAll, calAge_yr] = multiMatcalQ(age_kyr, ageErr_kyr, 1:length(age_kyr), S);

% Randomly choose which scenario to use for each run
numRuns        = S.numruns;
runScenarioIdx = randsample(1:length(scenarios), numRuns, 'true');

% Create S.numruns possible sedimentation histories, by random sampling
for iRun = 1:numRuns
    % Find which dates are in the random scenario
    inScenarioLog   = ismember(string(label), scenarios{runScenarioIdx(iRun)});
    scenarioDateIdx = find(inScenarioLog == 1);

    % Get ages' probability vectors
    ageProb = ageProbAll(:,scenarioDateIdx);

    % Get depths in this scenario
    scenario_depth_cm = depth_cm(scenarioDateIdx);

    % Assign a single age to each date, either using the modes or by random
    % sampling from the probability vectors
    if S.useModes
        % --- Using modes
        sampledAge_yr = ageModes{runScenarioIdx(iRun)};
        sampledAge_yr = sampledAge_yr';
    else
        % --- Using random sampling

        % Initialise
        sampledAge_yr = NaN(1,length(scenarioDateIdx));

        % Run through each date and randomly choose a year from its
        % probability vector, ensuring that each age is older than the age
        % directly shallower than it
        i = 1;

        % Set up how to run through each date (ascending, descending, random)
        shallow2deepOrder = 1:length(scenarioDateIdx);
        deep2shallowOrder = fliplr(shallow2deepOrder);
        randomOrder       = randperm(length(scenarioDateIdx));

        % Initialise the rejection counter
        numRejected = 0;

        if sampleOrder == "rejection"
            % --- Joint rejection sampling (reference sampler)
            % Draw every date independently from its full calibrated PDF and
            % keep the run only if the ages strictly increase with depth.
            % Drawing is vectorised: each date's PDF is turned into a
            % cumulative distribution once, then S.rejectionBatchSize candidate
            % runs are drawn at a time by mapping uniform random numbers through
            % it (inverse-CDF sampling via discretize)

            % Count the dates in this scenario
            numDates = length(scenarioDateIdx);

            % Turn each date's PDF (above S.pdfMinVal) into the edges of its CDF
            cdfEdges    = cell(1, numDates);
            candAges_yr = cell(1, numDates);
            for j = 1:numDates
                aboveMinProbLog  = ageProb(:,j) >= S.pdfMinVal;
                keptAgeProb      = ageProb(aboveMinProbLog,j);
                keptCdf          = cumsum(keptAgeProb)./sum(keptAgeProb);
                cdfEdges{j}      = [0; keptCdf];
                cdfEdges{j}(end) = 1;   % Guard against cumsum/sum rounding leaving it below 1
                candAges_yr{j}   = calAge_yr(aboveMinProbLog);
            end

            % Draw batches of candidate runs until one increases with depth
            foundMonotone = false;
            numBatches    = 0;
            while ~foundMonotone
                numBatches = numBatches + 1;
                if numBatches > S.rejectionMaxBatches
                    % Rejection sampling is not feasible for this core (the
                    % calibrated PDFs overlap too much), so treat the core as
                    % rejected rather than stopping the whole pipeline
                    warning("oneCoreRSR: no monotone run found for core " + corename + ...
                        " after " + num2str(S.rejectionBatchSize.*S.rejectionMaxBatches) + " candidates - " + ...
                        "rejection sampling is not feasible, core skipped")
                    nSRcounts = [];
                    agediffs  = [];
                    return
                end

                % Draw a batch of candidate runs by mapping uniform random
                % numbers through each date's CDF
                candRuns_yr  = NaN(S.rejectionBatchSize, numDates);
                uniformDraws = rand(S.rejectionBatchSize, numDates);
                for j = 1:numDates
                    cdfBinIdx        = discretize(uniformDraws(:,j), cdfEdges{j});
                    candRuns_yr(:,j) = candAges_yr{j}(cdfBinIdx);
                end

                % Find the first candidate run whose ages increase with depth
                monotoneLog      = all(diff(candRuns_yr, 1, 2) > 0, 2);
                firstMonotoneIdx = find(monotoneLog, 1);
                if ~isempty(firstMonotoneIdx)
                    sampledAge_yr = candRuns_yr(firstMonotoneIdx,:);
                    foundMonotone = true;
                end
            end

            %%% RESTRICTION CODE (applied after the draw, shallow to deep):
            % Drop any date within minAgeDiff of the last kept date
            if minAgeDiff > 0
                lastKeptAge_yr = sampledAge_yr(1);
                for j = 2:numDates
                    if sampledAge_yr(j) - lastKeptAge_yr < minAgeDiff
                        sampledAge_yr(j) = NaN;
                    else
                        lastKeptAge_yr = sampledAge_yr(j);
                    end
                end
            end
        else
            while i <= length(scenarioDateIdx)
                % Choose which date to sample at this step, following sampleOrder
                if sampleOrder == "shallow2deep"
                    iDate = shallow2deepOrder(i);
                elseif sampleOrder == "deep2shallow"
                    iDate = deep2shallowOrder(i);
                elseif sampleOrder == "random"
                    iDate = randomOrder(i);
                end

                % Set up age probabilities so that monotonicity is forced
                if i == 1
                    % Sample the first date from its PDF above S.pdfMinVal
                    aboveMinProbLog      = ageProb(:,iDate) >= S.pdfMinVal;
                    possAges_yr          = calAge_yr(aboveMinProbLog);
                    possAgeProbs         = ageProb(aboveMinProbLog,iDate);
                    sampledAge_yr(iDate) = randsample(possAges_yr, 1, true, possAgeProbs);
                else
                    % Find the oldest age already sampled for a shallower date
                    % (this date must be older than it)
                    minAge_yr = max(sampledAge_yr(1:iDate-1));
                    if isempty(minAge_yr) || isnan(minAge_yr)
                        minAge_yr = -Inf;
                    end

                    % Find the youngest age already sampled for a deeper date
                    % (this date must be younger than it)
                    maxAge_yr = min(sampledAge_yr(iDate+1:end));
                    if isempty(maxAge_yr) || isnan(maxAge_yr)
                        maxAge_yr = Inf;
                    end

                    % Get the calibrated ages between those two limits and
                    % their probabilities
                    possibleAgesLog = calAge_yr > minAge_yr & calAge_yr < maxAge_yr;
                    possAges_yr     = calAge_yr(possibleAgesLog);
                    possAgeProbs    = ageProb(possibleAgesLog,iDate);

                    % Normalise the probabilities and drop those below S.pdfMinVal
                    possAgeProbs    = possAgeProbs./sum(possAgeProbs);
                    aboveMinProbLog = possAgeProbs >= S.pdfMinVal;
                    possAges_yr     = possAges_yr(aboveMinProbLog);
                    possAgeProbs    = possAgeProbs(aboveMinProbLog);

                    % Choose the age for this date
                    if isempty(possAges_yr) && (isinf(minAge_yr) || isinf(maxAge_yr))
                        % Window runs off the end of the calibration axis: reset
                        % the run (original intent)
                        i = 1;
                        sampledAge_yr = NaN(1,length(scenarioDateIdx));
                        continue
                    elseif isempty(possAges_yr)
                        % Squeezed between two accepted neighbours: no room for this date
                        sampledAge_yr(iDate) = NaN;
                        numRejected = numRejected + 1;
                    elseif isscalar(possAges_yr)
                        sampledAge_yr(iDate) = possAges_yr;
                    else
                        sampledAge_yr(iDate) = randsample(possAges_yr, 1, true, possAgeProbs);
                    end

                    %%% RESTRICTION CODE
                    % If next age is not more than x years greater than
                    % previous accepted age, reject it and move to sample from
                    % next radiocarbon date (x stored in input minAgeDiff)
                    if sampledAge_yr(iDate) - minAge_yr < minAgeDiff || ...
                            maxAge_yr - sampledAge_yr(iDate) < minAgeDiff
                        sampledAge_yr(iDate) = NaN;
                    end
                end
                i = i+1;
            end
        end
    end

    % Find out which ages were used and which weren't
    usedLog      = ~isnan(sampledAge_yr);
    run_depth_cm = scenario_depth_cm(usedLog);
    run_age_yr   = sampledAge_yr(usedLog);
    %datesUsedStore(runN,:) = datesUsed;

    % Calculate the average SR of this run in cm/kyr (depths in cm, ages in yr)
    run_aveSR = (run_depth_cm(end)-run_depth_cm(1))./(run_age_yr(end)-run_age_yr(1)).*1000;

    % --- Calculate the SR of each pair of ages, in cm/kyr
    run_dt_yr = diff(run_age_yr);
    run_dx_cm = diff(run_depth_cm)';
    run_SRs   = run_dx_cm./run_dt_yr.*1000;

    % Normalise the SRs by this run's average SR or by the scenario's mean SR
    if S.normWithRunAve
        normSRs = run_SRs./run_aveSR;
    else
        normSRs = run_SRs./scenario_meanSR(runScenarioIdx(iRun));
    end

    % Find normalising value based on number of runs
    weightingNormaliser = numRuns;

    % Add normSRs to vector to count them (with their weighting)
    nSRinfo = [normSRs; run_dx_cm; run_dt_yr];   % nSR, dx in cm, dt in yr

    % Store all nSR info in a single array, with NaNs separating info from
    % different runs (instead of NaN in age-diff row, the first age used
    % is stored, so that the whole time history can be created from this)
    if iRun == 1
        nSRcounts = [[NaN; min(run_depth_cm); min(run_age_yr)], nSRinfo];
        agediffs  = run_dt_yr;
    else
        nSRcounts = cat(2, nSRcounts, [[NaN; min(run_depth_cm); min(run_age_yr)], nSRinfo]);
        agediffs  = cat(2, agediffs, run_dt_yr);
    end
end
end
