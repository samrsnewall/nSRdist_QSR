function [invSRvals, invSRprobs, aveSR, reversalPairs, numDatePairs, ageModes_yr, ...
    sedLength, ageLength_yr, MSI_byage, MSI_bydepth, IDpairs, agediffV] = ...
    scenariopdfNorm(depth_cm, scenarioDateIdx, label, ageProb, calAge_yr, IDpairs, ...
    agediffV, S, plotfigs)
% This function takes the dates of one age-depth scenario (scenarioDateIdx)
% and gets the PDF of the age difference of each pair of neighbouring dates
% with agediffcalc. A pair is a reversal if the probability of an age
% difference of zero or less is S.reversalCriteria or more; if any pair is,
% the function returns straight away with invSRvals, invSRprobs and aveSR
% set to NaN (the other outputs are filled in), so that
% scenariosDealWithReversals can deal with it. Otherwise it works out the
% scenario's average SR and, only if S.pdfMethod, its PDF of normalised
% inverse SR, summed over the pairs (weighted by their depth differences
% if S.weighting is "depth").
%
% depth_cm and label hold every date in the core, and ageProb the
% calibrated age PDFs of this scenario's dates only (one column per date,
% on the calAge_yr grid). Depths and sedLength (the depth span) in cm;
% calibrated ages, ageModes_yr, ageLength_yr (the span of the mode ages)
% and age differences in yr; aveSR in cm/kyr (the depth span over the
% mode-age span in kyr). MSI_byage is ageLength_yr divided by the number of
% dates (yr per date) and MSI_bydepth is sedLength divided by the number of
% dates (cm per date). Each pair's inverse SR is in yr/cm and is normalised
% by multiplying by aveSR/1000, so invSRvals is unit-free; invSRvals and
% invSRprobs are empty unless S.pdfMethod. reversalPairs holds 1 for each
% pair that is a reversal and 0 otherwise. IDpairs and agediffV are the
% caches of date pairs (by LabID) and their age-difference PDFs, passed to
% agediffcalc and returned with any new pairs added; plotfigs turns on the
% diagnostic plots.
%
% S fields used: reversalCriteria, pdfMethod, weighting (only if
% S.pdfMethod), and pdfMinVal through agediffcalc

%% Convert the calibrated age grid to kyr
% Used for the average SR below
calAge_kyr = calAge_yr ./ 1000;

%% Find the mode ages and the spans of the scenario
% Take the mode age of each date (the data are already sorted by depth)
ageModes_yr = zeros(length(scenarioDateIdx), 1);
for i = 1:length(scenarioDateIdx)
    [~, modeIdx]   = max(ageProb(:,i));
    ageModes_yr(i) = calAge_yr(modeIdx);
end

% Get the depths of the dates in this scenario, the age span between the
% mode ages of the first and last dates, and the depth span
scenario_depth_cm = depth_cm(scenarioDateIdx);
ageLength_yr      = ageModes_yr(end) - ageModes_yr(1);
sedLength         = scenario_depth_cm(end) - scenario_depth_cm(1);

% Calculate the mean sampling interval (MSI) by age (yr per date) and by
% depth (cm per date), dividing by the number of dates
MSI_byage   = ageLength_yr / length(scenarioDateIdx);
MSI_bydepth = sedLength / length(scenarioDateIdx);

%% Get the age-difference PDF of each pair of neighbouring dates
% agediffcalc takes a pair from the cache (IDpairs, agediffV) if it has
% been worked out before, and adds it to the cache if not
numDatePairs   = length(scenarioDateIdx) - 1;
scenario_label = label(scenarioDateIdx);
[dtVals_yr, dtProbs, IDpairs, agediffV] = agediffcalc(ageProb, calAge_yr, numDatePairs, ...
    scenario_label, IDpairs, agediffV, S);

% Calculate the depth difference of each pair
dx_cm = diff(scenario_depth_cm);

%% Find which pairs are reversals
% Call a pair a reversal if the probability of an age difference of zero
% or less is S.reversalCriteria or more (an arbitrary threshold)
reversalPairs = zeros(1, numDatePairs);
for iPair = 1:numDatePairs
    nonPositiveDtLog = dtVals_yr{iPair} <= 0;
    reversalProb     = sum(dtProbs{iPair}(nonPositiveDtLog));
    if reversalProb >= S.reversalCriteria
        reversalPairs(iPair) = 1;
    end
end

% Stop here if any pair is a reversal, so that scenariosDealWithReversals
% can remove a date from this scenario or split it
if sum(reversalPairs) > 0
    invSRvals  = NaN;
    invSRprobs = NaN;
    aveSR      = NaN;
    return
end

%% Keep only the positive age differences
% Drop any probability at age differences of zero or less, and renormalise
% what is left of each pair's PDF to sum to 1
dtProbsPos   = cell(1, numDatePairs);
dtValsPos_yr = cell(1, numDatePairs);
for iPair = 1:numDatePairs
    positiveDtLog       = dtVals_yr{iPair} > 0;
    dtValsPos_yr{iPair} = dtVals_yr{iPair}(positiveDtLog);
    dtProbsPos{iPair}   = dtProbs{iPair}(positiveDtLog);
    dtProbsPos{iPair}   = dtProbsPos{iPair} ./ sum(dtProbsPos{iPair}, 'all');
end

%% Convert the age differences to inverse SR
% Divide each pair's age differences by its depth difference (so the
% inverse SR here is in yr/cm, not kyr/cm), and renormalise its PDF to
% unit area
invSRvals_byPair  = cell(1, numDatePairs);
invSRprobs_byPair = cell(1, numDatePairs);
invSRAUC          = NaN(1, numDatePairs);
for iPair = 1:numDatePairs
    invSRvals_byPair{iPair}  = dtValsPos_yr{iPair} ./ dx_cm(iPair);
    invSRAUC(iPair)          = trapz(invSRvals_byPair{iPair}, dtProbsPos{iPair});
    invSRprobs_byPair{iPair} = dtProbsPos{iPair} ./ invSRAUC(iPair);
end

%% Calculate the average SR of the scenario
% Divide the depth span (cm) by the span of the mode ages of the first and
% last dates (kyr), giving cm/kyr. These are the same spans as sedLength
% and ageLength_yr above, worked out again here
[~, topModeIdx]    = max(ageProb(:,1));
[~, bottomModeIdx] = max(ageProb(:,end));
topAgeMode_kyr     = calAge_kyr(topModeIdx);
bottomAgeMode_kyr  = calAge_kyr(bottomModeIdx);
ageLength_kyr      = bottomAgeMode_kyr - topAgeMode_kyr;
sedLengthForSR     = scenario_depth_cm(end) - scenario_depth_cm(1);
aveSR              = sedLengthForSR ./ ageLength_kyr;

%% Normalise the inverse SR by the average SR
% Multiply each pair's inverse SR (yr/cm) by aveSR, the same as dividing it
% by the inverse of aveSR, with *(1/1000) to put aveSR in cm/yr, so the
% normalised inverse SR is unit-free
normInvSRvals_byPair = cellfun(@(x) x * (aveSR .* (1/1000)), invSRvals_byPair, 'un', 0);

%% Plot each pair's PDF
% Diagnostic plots of every pair's PDF of inverse SR (in kyr/cm) and of
% normalised inverse SR
if plotfigs == 1
    pdfCreation = figure("Name", "pdfCreation");

    subplot(4,1,1)
    for iPair = 1:numDatePairs
        plot(invSRvals_byPair{1,iPair} ./ 1000, invSRprobs_byPair{1,iPair}, 'LineWidth', 1)
        hold on
    end
    ylabel("Probability")
    xlabel("Inverse of Sed Rate ky/cm")

    subplot(4,1,2)
    hold on
    for iPair = 1:numDatePairs
        plot(normInvSRvals_byPair{1,iPair}, invSRprobs_byPair{1,iPair})
    end
    ylabel("Probability")
    xlabel("Normalised Inverse Sed Rate")
end

%% Build the scenario's PDF of normalised inverse SR
% Only for the analytical PDF method (S.pdfMethod): put every pair's PDF on
% one grid, weight the pairs, and sum them
if S.pdfMethod

    % --- Interpolate each pair's PDF onto a common grid
    % Set up a grid of normalised inverse SR with a constant spacing that
    % runs through multiples of the spacing (e.g. 15.0, 15.1, 15.2 instead
    % of 15.032, 15.132, 15.232), so the pairs' PDFs can be summed
    spacing          = 0.00005;
    normInvSRvalsAll = vertcat(normInvSRvals_byPair{:});
    lowerBound       = floor(min(normInvSRvalsAll * (1/spacing))) * spacing;
    upperBound       = max(normInvSRvalsAll);
    invSRvals        = lowerBound:spacing:(upperBound + spacing * 0.9);

    % Interpolate each pair's PDF onto the part of the grid within its
    % range, and normalise the interpolated PDF to sum to 1
    interpInvSRvals_byPair  = cell(1, numDatePairs);
    interpInvSRprobs_byPair = cell(1, numDatePairs);
    for iPair = 1:numDatePairs
        if ~isempty(normInvSRvals_byPair{iPair})
            % Find where on the grid this pair's values fit
            minNormInvSR = min(normInvSRvals_byPair{iPair});
            maxNormInvSR = max(normInvSRvals_byPair{iPair});
            firstGridIdx = find(invSRvals > minNormInvSR, 1, 'first');
            lastGridIdx  = find(invSRvals < maxNormInvSR, 1, 'last');

            % Interpolate within this part of the grid, and normalise
            interpInvSRvals_byPair{iPair}  = invSRvals(firstGridIdx:lastGridIdx);
            interpProbs                    = interp1(normInvSRvals_byPair{iPair}, ...
                invSRprobs_byPair{iPair}, invSRvals(firstGridIdx:lastGridIdx));
            interpInvSRprobs_byPair{iPair} = interpProbs ./ sum(interpProbs);
        end
    end

    % --- Line up the pairs' PDFs on the grid and weight them
    % Put each pair's interpolated PDF in its own column of a matrix on the
    % full grid, starting at the grid value that matches its first value (to
    % within 0.00001)
    invSRprobsArray = zeros(length(invSRvals), numDatePairs);
    for iPair = 1:numDatePairs
        if ~isempty(interpInvSRvals_byPair{iPair})
            startIdx    = find(abs(invSRvals - interpInvSRvals_byPair{iPair}(1)) <= 0.00001);
            numGridVals = length(interpInvSRvals_byPair{iPair}(:));
            invSRprobsArray(startIdx:startIdx + numGridVals - 1, iPair) = ...
                interpInvSRprobs_byPair{iPair}(:);
        end
    end

    % Weight each pair's PDF by its depth difference if S.weighting is
    % "depth"
    % FIXME: calcData does not set S.weighting (weighting is fitS.weighting
    % in fitData), so with S.pdfMethod true the line below stops with an
    % error (as noted in oneCoreScenarios)
    if S.weighting == "depth"
        weightedInvSRprobsArray = invSRprobsArray .* dx_cm';
    else
        % "age" weighting is not yet implemented (this cannot weight by
        % age), and "none" uses equal weights
        weightedInvSRprobsArray = invSRprobsArray;
    end

    % Plot each pair's weighted PDF on the grid
    if plotfigs == 1
        figure(pdfCreation)
        subplot(4,1,3)
        hold on
        plot(invSRvals, weightedInvSRprobsArray)
        xlabel("Individual invSR Ratios weighted")
        ylabel("Probability")
    end

    % --- Sum the pairs' PDFs
    % Sum the weighted PDFs of all the pairs, and normalise to unit area
    summedInvSRprobs = sum(weightedInvSRprobsArray, 2);
    summedAUC        = trapz(invSRvals, summedInvSRprobs);
    invSRprobs       = summedInvSRprobs ./ summedAUC;

    % Plot the scenario's PDF
    if plotfigs == 1
        figure(pdfCreation)
        subplot(4,1,4)
        hold on
        plot(invSRvals, invSRprobs)
        ylabel("Probability")
        xlabel("Inverse of Sed Rate Ratio (Weighted by depth)")
    end

else
    % No PDF without S.pdfMethod
    invSRvals  = [];
    invSRprobs = [];
end

end
