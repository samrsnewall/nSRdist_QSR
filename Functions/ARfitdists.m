function[outS] = ARfitdists(dataTCol, x, chooseLog, weightRepDataDP, weightRepInflator, countDivisor, applyLin2014AgeFilter, fitS)
% This function pools the nSR data of the chosen cores (chooseLog) for one
% method, weights them by depth, by age or not at all (fitS.weighting), and
% fits a 2-component mixed log-normal (MLN), a log-normal (LN), a gamma (Gam)
% and an inverse gamma (invGam) to the pooled data, giving each fitted pdf
% on the grid x and transformed to log(nSR). If applyLin2014AgeFilter is
% true it first keeps only the pairs whose dt_yr is at least the lower value
% of fitS.Lin2014AgeFilter and below the upper value (fitData does this for
% the <method>Lin fits). nSR is unit-free (SR is in cm/kyr if
% fitS.non_normalized_SR); dx in cm; dt in yr.
%
% dataTCol is one dataT column (one 3-row nSR matrix per core: nSR, dx_cm,
% dt_yr, each run starting with a NaN column). weightRepDataDP and
% weightRepInflator set how the weighted data are replicated (see
% makeWeightedReplicates). countDivisor is the number of runs per core the
% method is taken to have: each pair is also weighted by countDivisor
% divided by its own core's number of runs, so every core counts as
% countDivisor runs' worth whatever number it has (see "Weight each pair by
% its core's number of runs" below), and outS.numCpairs is the number of
% nSR values per run summed over the cores, which the fitters use as the
% effective number of observations. outS holds the MLN, LN, Gam and
% invGam results (each with .fitInfo, the .nSR and .lnSR pdfs with their
% mean and variance, and .chiStats if fitS.run_chi2gof), the data the
% fitters saw (.weightedC), .numCpairs, .Lin2014AgeFiltering, the weighted
% bin counts of dt_yr (.agediffsWC), the summed dx (cm) and dt (yr) of the
% chosen cores per run (.sedLength, .sedTimeSpan), and .TM.
%
% fitS fields used: non_normalized_SR, Lin2014AgeFilter,
% weighting, weightingPipeline, resampleData, maxFitPoints, mlnReps,
% run_chi2gof, calculateTM (and those chi2_dataVStwopdfVECs and
% TMcalculationAR read)

%% Combine all counts into one array
% Put the nSR matrices of the chosen cores side by side
nSRcountsArray         = countsCell2Array(dataTCol, chooseLog);
nSRcountsArrayAsLoaded = nSRcountsArray;   % Copy before the checks below (not used further on)

% Count the runs of each chosen core with data (its NaN run-start columns),
% in the order countsCell2Array put the cores side by side, so that each
% pair can be weighted by its core's number of runs below
runsPerCore = [];
for i = 1:length(dataTCol)
    if chooseLog(i) == 1 && ~isempty(dataTCol{i})
        runsPerCore = [runsPerCore, sum(isnan(dataTCol{i}(1,:)))]; %#ok<AGROW>
    end
end

% Remove negative SRs, keeping the NaN run-start columns (ideally there
% shouldn't be any)
if sum(nSRcountsArray(1,:) < 0) ~= 0
    warning("There are negative sed rates!... being removed")
    notNegativeLog = nSRcountsArray(1,:) >= 0 | isnan(nSRcountsArray(1,:));
    nSRcountsArray = nSRcountsArray(:, notNegativeLog);
end

% Remove SRs of 0, keeping the NaN run-start columns (ideally there
% shouldn't be any)
if sum(nSRcountsArray(1,:) == 0) ~= 0
    warning("There are SRs of 0!... being removed")
    positiveLog    = nSRcountsArray(1,:) > 0 | isnan(nSRcountsArray(1,:));
    nSRcountsArray = nSRcountsArray(:, positiveLog);
end

% If wanted, use absolute SR (cm/kyr) instead of nSR. Row 2 is dx in cm
% and row 3 is dt in yr, so the ratio is multiplied by 1000 to give cm/kyr
if fitS.non_normalized_SR
    dx_cm = nSRcountsArray(2,:);
    dt_yr = nSRcountsArray(3,:);
    SRs   = dx_cm./dt_yr.*1000;

    % Keep the run-start columns as NaN, and put the SRs in row 1
    runStartLog         = isnan(nSRcountsArray(1,:));
    SRs(runStartLog)    = NaN;
    nSRcountsArray(1,:) = SRs;
end

% Keep a copy before the Lin2014 age filter, for the TM
nSRcountsArrayNoAgeFilter = nSRcountsArray;

% Apply the age filter as done by Lin2014 if desired: keep the pairs whose
% dt_yr is at least the lower value of fitS.Lin2014AgeFilter and below the
% upper value (the lower bound is inclusive so that BMedian pairs with a dt
% of exactly 500 yr, common with whole-year median ages, are kept), and the
% NaN run-start columns
if applyLin2014AgeFilter
    maxDt_yr       = max(fitS.Lin2014AgeFilter);
    minDt_yr       = min(fitS.Lin2014AgeFilter);
    inAgeFilterLog = nSRcountsArray(3,:) < maxDt_yr & nSRcountsArray(3,:) >= minDt_yr;
    runStartLog    = isnan(nSRcountsArray(1,:));
    lin2014KeepLog = inAgeFilterLog | runStartLog;
    nSRcountsArray = nSRcountsArray(:,lin2014KeepLog);
end

%% Weight each pair by its core's number of runs
% Pooling every run of every core weights each core by its number of runs,
% and those differ between cores for BSamp (runs with too few ages are
% skipped). So each pair gets a run weight of countDivisor divided by its
% core's number of runs: one for a core with countDivisor runs, so a method
% whose cores all have countDivisor runs (RSRx, rejection, BMedian) is
% unchanged, and every core counts as countDivisor runs' worth whatever
% number it has. The core of each column is found from the run-start
% columns, which every filter above keeps (one per run), in the order the
% cores were put side by side
runStartLog = isnan(nSRcountsArray(1,:));
runIdx      = cumsum(runStartLog);                              % Which run each column is in
runCoreIdx  = discretize(runIdx, [0, cumsum(runsPerCore)] + 0.5); % Which core each run is in
runWeights  = countDivisor./runsPerCore(runCoreIdx);

%% Apply weighting
% Convert the weighted nSR counts to a single dimension array of counts
% where the number of counts is representative of their weighting

% Set up nSR counts and weighting vectors without the NaN run-start columns
% (nSRs are SRs in cm/kyr if fitS.non_normalized_SR)
notNanLog      = ~isnan(nSRcountsArray(1,:));
nSRs           = nSRcountsArray(1,notNanLog);
dx_cm          = nSRcountsArray(2,notNanLog);
dt_yr          = nSRcountsArray(3,notNanLog);
pairRunWeights = runWeights(notNanLog);

% Set up the weights according to the weighting choice, times the run weight
if fitS.weighting == "none"
    weights = pairRunWeights;
elseif fitS.weighting == "depth"
    weights = dx_cm.*pairRunWeights;   % Weight by depth difference
elseif fitS.weighting == "age"
    weights = dt_yr.*pairRunWeights;   % Weight by time span
else
    error("Weighting type not recognized, must be 'none', 'depth', or 'age'")
end

% Record the number of nSR calculations per run, summed over the cores
% (used by the fitters as the effective number of observations), and the
% per-method age-filter choice, since fitS is saved once for all methods
outS.numCpairs           = sum(pairRunWeights)./countDivisor;
outS.Lin2014AgeFiltering = applyLin2014AgeFilter;

% Turn the weighted data into the vector the fitters see (weightedC). Two
% pipelines that draw from the same distribution; fitS.weightingPipeline
% chooses:
%   "replicate" - build the full replicated dataset (every value repeated
%                 according to its weight), then, if it is bigger than
%                 fitS.maxFitPoints, take a random subsample of that size
%   "sample"    - draw the fitS.maxFitPoints subsample directly from the
%                 replicate counts without building the replicated dataset
%                 (same distribution, much faster on large datasets)
if fitS.weightingPipeline == "replicate"

    % Weight by replicating data according to weighting. With no weighting
    % and every core at countDivisor runs the data are used as they are;
    % with no weighting but unequal runs the run weights are replicated
    % like any other weights (rounded to whole replicates with
    % weightRepInflator, which is 1 for no weighting, so fractional run
    % weights are coarse unless it is raised)
    if fitS.weighting == "none" && all(pairRunWeights == 1)
        weightedC = nSRs;
    else
        weightedC = makeWeightedReplicates(nSRs, weights, weightRepDataDP, weightRepInflator);
    end

    % Resample data back to the original number of observations
    % (outS.numCpairs, the number per run) if desired
    if fitS.resampleData
        rng(3)
        allWeightedC = weightedC;
        weightedC    = randsample(allWeightedC, floor(outS.numCpairs), false);
    end

    % Cap the number of replicated points passed to the fitters: with many
    % cores and runs the replicated data reaches 1e8+ points, and fitgmdist's
    % EM iterations over that take hours per fit. A random subsample of a few
    % million points gives the same fits, and BICeff is unaffected in
    % expectation because the NLL and its replication divisor shrink by the
    % same factor
    if length(weightedC) > fitS.maxFitPoints
        rng(4)
        weightedC = randsample(weightedC, fitS.maxFitPoints, false);
    end

elseif fitS.weightingPipeline == "sample"

    % Draw the capped subsample straight from the weighted data (the same
    % draw the "replicate" pipeline makes, without writing out the replicates)
    if fitS.resampleData
        error("fitS.resampleData is only implemented in the ""replicate"" weighting pipeline")
    end
    rng(4)
    if fitS.weighting == "none" && all(pairRunWeights == 1)
        weightedC = nSRs;
        if length(weightedC) > fitS.maxFitPoints
            weightedC = randsample(weightedC, fitS.maxFitPoints, false);
        end
    else
        weightedC = sampleWeightedReplicates(nSRs, weights, weightRepDataDP, weightRepInflator, ...
            fitS.maxFitPoints);
    end

else
    error("fitS.weightingPipeline must be ""replicate"" or ""sample""")
end

%% Fit the mixed log-normal (MLN) to these counts
% Fit a 2-component MLN and get its pdf on the grid x
[mixLogNormPDF, ~, MLN.fitInfo] = fitMixLogNorm(weightedC, x, 2, fitS.mlnReps, outS.numCpairs);
MLN.nSR.x  = mixLogNormPDF(:,1);
MLN.nSR.px = mixLogNormPDF(:,2);

%% Fit the log-normal (LN) to these counts
% Fit a 1-component MLN (a log-normal) and get its pdf on the grid x
[logNormPDF, ~, LN.fitInfo] = fitMixLogNorm(weightedC, x, 1, fitS.mlnReps, outS.numCpairs);
LN.nSR.x  = logNormPDF(:,1);
LN.nSR.px = logNormPDF(:,2);

%% Fit the gamma (Gam) to these counts
% Estimate the gamma fit parameters and get its pdf on the grid x
[Gam.nSR.px, Gam.fitInfo] = fitGamma(weightedC, x, outS.numCpairs);
Gam.nSR.x = x';

%% Fit the inverse gamma (invGam) to these counts
% Fit a gamma to the inverse of the counts and get the inverse gamma pdf on
% the grid x
[invGam.nSR.px, ~, invGam.fitInfo] = fitInvGamma(weightedC, x, outS.numCpairs);
invGam.nSR.x = x';

%% Save other useful information
% Save the data the fitters saw
outS.weightedC = weightedC;

% Calculate the weighted bin counts of dt_yr (bin edges in yr)
dtBinEdges_yr   = 0:100:10000;
outS.agediffsWC = makeWeightedBinCounts(dt_yr, weights, dtBinEdges_yr);

% Calculate the summed dx (cm) and dt (yr) of the chosen cores per run:
% each core's total divided by its own number of runs, summed over the
% cores (the run weights do this)
outS.sedLength   = sum(dx_cm.*pairRunWeights, "omitnan")./countDivisor;
outS.sedTimeSpan = sum(dt_yr.*pairRunWeights, "omitnan")./countDivisor;

%% Calculate the mean and variance of each pdf
[invGam.nSR.mu, invGam.nSR.var] = muVarPDFVec(invGam.nSR);
[Gam.nSR.mu, Gam.nSR.var]       = muVarPDFVec(Gam.nSR);
[LN.nSR.mu, LN.nSR.var]         = muVarPDFVec(LN.nSR);
[MLN.nSR.mu, MLN.nSR.var]       = muVarPDFVec(MLN.nSR);

%% Transform the pdfs to log(nSR) space
% Transform each pdf to log(nSR), calculate its mean and variance there,
% and record its number of parameters and name (read by
% chi2_dataVStwopdfVECs and the results scripts)
logFH = @(x) log(x);

% --- MLN
[MLN.lnSR.x, MLN.lnSR.px]   = px_to_pfx(MLN.nSR.x, MLN.nSR.px, logFH);
[MLN.lnSR.mu, MLN.lnSR.var] = muVarPDFVec(MLN.lnSR);
MLN.lnSR.numParams          = 5;
MLN.lnSR.pdfName            = "2 Component Mix Log Normal";

% --- invGam
[invGam.lnSR.x, invGam.lnSR.px]   = px_to_pfx(invGam.nSR.x, invGam.nSR.px, logFH);
[invGam.lnSR.mu, invGam.lnSR.var] = muVarPDFVec(invGam.lnSR);
invGam.lnSR.numParams             = 2;
invGam.lnSR.pdfName               = "Inverse Gamma";

% --- Gam
[Gam.lnSR.x, Gam.lnSR.px]   = px_to_pfx(Gam.nSR.x, Gam.nSR.px, logFH);
[Gam.lnSR.mu, Gam.lnSR.var] = muVarPDFVec(Gam.lnSR);
Gam.lnSR.numParams          = 2;
Gam.lnSR.pdfName            = "Gamma";

% --- LN
[LN.lnSR.x, LN.lnSR.px]   = px_to_pfx(LN.nSR.x, LN.nSR.px, logFH);
[LN.lnSR.mu, LN.lnSR.var] = muVarPDFVec(LN.lnSR);
LN.lnSR.numParams         = 2;
LN.lnSR.pdfName           = "LogNorm";

%% Run the chi-squared goodness-of-fit tests if desired
if fitS.run_chi2gof
    % Make a column cell-vector that holds all pdfs to test
    pdfs = {LN.lnSR; MLN.lnSR; Gam.lnSR; invGam.lnSR};

    % Test log(weightedC) against each pdf, starting from 20 bins
    % TODO: fitS.dispChi2 is forced false here, so fitData's setting has no
    % effect on these tests, and the 20 initial bins could come from fitS
    % (fitS.chi2binN is not used here)
    fitS.dispChi2   = false;
    [h, p, chiStat] = chi2_dataVStwopdfVECs(log(weightedC), outS.numCpairs, 20, pdfs, fitS);

    % Store the chiStats in the data structures
    for i = 1:size(pdfs,1)
        chiStat{i}.h = h(i);
        chiStat{i}.p = p(i);
    end
    LN.chiStats     = chiStat{1};
    MLN.chiStats    = chiStat{2};
    Gam.chiStats    = chiStat{3};
    invGam.chiStats = chiStat{4};
end

%% Calculate TM
% Calculate the transition matrix from the counts before the Lin2014 age
% filter, if desired. Not sure from Lin2014 whether the TM was calculated
% with or without the nSR counts from outside the dt range. If without, how
% does it deal with the gaps?
% TODO: TMcalculationAR does not use the run weights, so its counts weight
% each core by its number of runs
if fitS.calculateTM
    [~,~,TM,~] = TMcalculationAR(nSRcountsArrayNoAgeFilter, fitS);
else
    TM = [];
end

%% Put structures into output structure
outS.LN     = LN;
outS.MLN    = MLN;
outS.Gam    = Gam;
outS.invGam = invGam;
outS.TM     = TM;

end
