%% fitData - stage 2 of the nSRdist pipeline
% This script loads a calcData output file and fits four candidate
% probability distributions to the nSR data from each nSR estimation method
% (see below), and optionally individual-run variants. It saves all the fit
% results for downstream plotting and analysis. So calcData saves
% Results/dataT_*.mat, and fitData then loads that file, fits distributions
% and saves d_*.mat
%
% The distributions fitted by ARfitdists (pooled-run, or "AR", fits) are:
%   MLN    - 2-component Mixed Log-Normal, fitted by EM (fitgmdist)
%   LN     - Single-component Log-Normal (special case of MLN with 1
%            component)
%   Gam    - Gamma distribution fitted to nSR values
%   invGam - Inverse Gamma: a Gamma fitted to 1/nSR, then transformed back
% IRfitdists (individual-run, or "IR", fits, if fitS.fitIRs = true) fits the
% same four distributions to each of S.OneRun.numruns individual random
% samples (one run per core), giving a distribution of fit parameters whose
% spread reflects run-to-run variability
%
% The nSR estimation methods analysed (field names in d.S1) are:
%   <method>Lin         - Fitted with the Lin2014 age filter, for each
%                         method in fitS.lin2014Methods (BMedianLin,
%                         BSampLin)
%   RSR<dt><order>      - Sequential RSR, restriction applied while
%                         sampling, for every RSR column calcData made
%                         (RSR500 for the manuscript's shallow-to-deep
%                         order, RSR500d2s, RSR500random, ...)
%   <method>            - The draws of a method that samples every age
%                         first, fitted as drawn, for each in
%                         fitS.unthinnedMethods (rejection, BSamp, BMedian),
%                         no Lin2014 filter
%   <name>IR            - Individual-run fits of BSamp (BSampLinIR) and each
%                         RSR column (RSR500IR, ...), if fitS.fitIRs
%
% d.S1.chooseLog selects the cores that pass all of:
%   meanSR > fitS.minMeanSR  (excludes very low accumulation rate cores)
%   depths > fitS.minDepth   (excludes shallow-water cores)
%   |lat|  < fitS.maxLat     (excludes high-latitude cores)
%
% Saved to fileName2save:
%   d  - Struct containing all fit results, with fields:
%     d.dataFile - The calcData file that was fitted (load it for dataT)
%     d.S        - The calcData settings structure
%     d.x        - nSR grid (exp(d.logx)) on which PDFs are evaluated
%     d.logx     - Log-space grid for nSR
%     d.S1       - Results struct for the primary core subset:
%       d.S1.fitS       - Copy of the fitS settings used
%       d.S1.chooseLog  - Logical index of cores included in this subset
%       d.S1.numCores   - Number of cores in this subset
%       d.S1.<fit name> - ARfitdists (or IRfitdists) output for each fit
%                         listed under the nSR estimation methods above
%   BM - BIGMACS reference distribution (fitted MLN and histogram)
%
% The settings are in fitS, set in the "Set up fitting settings structure
% (fitS)" section; see the comment on each field there. A driver script can
% replace some of them, and the file to fit, by defining fitOverride before
% it runs this script; run standalone there is no fitOverride, so no effect
%
% It needs Functions/ (all helper functions, added to the path at startup)
% and Results/ (added to the path), which must contain the calcData output
% file

%% Add important paths
addpath('Functions')
addpath('Results')

%% Load the calcData results
filePath = "Results/dataT_QSR";

% Overrides from a driver script: if a calling script has defined
% fitOverride, its fields replace filePath (here) and fitS.weighting (below).
% No effect standalone
if exist('fitOverride', 'var')
    if isfield(fitOverride, 'filepath')
        filePath = fitOverride.filepath;
    end
end
load(filePath)

% FIXME: d and fitS are not cleared before they are filled, so running this
% script twice in one workspace carries fields from the first run (e.g. fits
% of methods no longer in the lists, or IR fits) into the second fit file

% Keep the file name rather than a copy of dataT, which is over 2 GB with
% every sampling order; scripts that need dataT load it from d.dataFile
d.dataFile = filePath;
d.S        = S;

% Get the file path without .mat, to build the fit file name from
filePathParts = split(filePath, ".mat");

%% Set up fitting settings structure (fitS)
% fitS controls all aspects of how distributions are fitted to the nSR data

% --- Core subset filters (applied to dataT)
fitS.minMeanSR = 8;    % Minimum mean SR (cm/kyr); cores below this are excluded
fitS.minDepth  = 1000; % Minimum water depth (m)
fitS.maxLat    = 40;   % Maximum absolute latitude (degrees)
% Note, filters (that may be different) are also applied in calcData

% --- Distribution fitting options
fitS.fitDists  = true;    % If true, fits all four distributions (MLN, LN, Gam, invGam)
fitS.weighting = "depth"; % How nSR estimates are weighted when fitting: "depth", "age", or "none"
if exist('fitOverride', 'var') && isfield(fitOverride, 'weighting')
    fitS.weighting = fitOverride.weighting;
end

% Set up the fit file name: calcData file name + "_fit_<weighting>weight"
fitSuffix     = "_fit_" + fitS.weighting + "weight";
fileName2save = filePathParts(1) + fitSuffix + ".mat";

% Save to a different folder if a driver script asks (fitOverride.saveFolder),
% keeping the same file name, so drivers fitting the same calcData file don't
% block each other
if exist('fitOverride', 'var') && isfield(fitOverride, 'saveFolder')
    [~, fitFileName] = fileparts(fileName2save);
    fileName2save    = fullfile(fitOverride.saveFolder, fitFileName + ".mat");
end

fitS.non_normalized_SR = false; % If true, fits absolute SR (cm/kyr) rather than nSR
fitS.mlnReps           = 3;     % Number of EM random restarts when fitting the Mixed Log-Normal;
                                %   more restarts reduce the chance of converging to a local
                                %   optimum

% --- Age filtering
fitS.Lin2014AgeFilter = [500 4000]; % [min max] age difference (yr) allowed between dated pairs
                                    %   when a method is fitted with applyLin2014AgeFilter = true;
                                    %   mirrors Lin2014 criterion

% --- Methods to fit (every RSR column calcData made is always fitted)
fitS.lin2014Methods   = ["BMedian", "BSamp"];
                                     % Fitted with the Lin2014 age filter, stored as <method>Lin
fitS.unthinnedMethods = ["rejection", "BSamp", "BMedian"];
                                     % Methods that sample every age first, fitted as drawn (no
                                     %   filter). BMedian is here for the unfiltered fit plotted
                                     %   behind BMedianLin in the manuscript Figure 2

% --- Transition Matrix
fitS.calculateTM = false;

% --- Pooled-run (AR) fitting parameters
fitS.invXbinEdges      = 0:0.1:15;    % Bin edges for inverse-SR histogram (used in invGam fitting)
fitS.resampleData      = false;       % If true, after weighting resamples data back to original N
                                      %   (removes inflation from replication)
fitS.maxFitPoints      = Inf;         % Cap on replicated points passed to the fitters (random
                                      %   subsample above this); Inf = no cap
fitS.weightingPipeline = "replicate"; % How the weighted data are prepared for the fitters:
                                      %   "replicate" builds the full replicated dataset then
                                      %   subsamples it to maxFitPoints (makeWeightedReplicates);
                                      %   "sample" draws the maxFitPoints subsample directly from
                                      %   the replicate counts (sampleWeightedReplicates) - same
                                      %   distribution, much faster on large datasets

% A driver script can change the cap and the pipeline; no effect standalone
if exist('fitOverride', 'var') && isfield(fitOverride, 'maxFitPoints')
    fitS.maxFitPoints = fitOverride.maxFitPoints;
end
if exist('fitOverride', 'var') && isfield(fitOverride, 'weightingPipeline')
    fitS.weightingPipeline = fitOverride.weightingPipeline;
end

% Set up the replication settings for the chosen weighting
if fitS.weighting == "depth"
    fitS.PooledRuns.weightRepDP       = 3;   % Decimal places to round data to in replication process
    fitS.PooledRuns.weightRepInflator = 0.5; % Multiplier on weights before replication (increase or
                                             %   reduce size of weights; 1 = no inflation)
elseif fitS.weighting == "age"
    fitS.PooledRuns.weightRepDP       = 3;
    fitS.PooledRuns.weightRepInflator = 0.001;
else
    fitS.PooledRuns.weightRepDP       = 3;
    fitS.PooledRuns.weightRepInflator = 1;
    fitS.OneRun.weightRepInflator     = 1;
end

% --- Individual-run (IR) fitting (BSampIR, RSR0IR, etc.)
fitS.fitIRs = false; % If true, fits distributions to each run separately (slower; gives
                     %   variability estimate)

% A driver script can switch IR fitting on or off; no effect standalone
if exist('fitOverride', 'var') && isfield(fitOverride, 'fitIRs')
    fitS.fitIRs = fitOverride.fitIRs;
end
if fitS.fitIRs
    fitS.OneRun.numruns = 1000; % Number of individual runs to sample and fit

    % A driver script can ask for fewer (the regression check does); no effect standalone
    if exist('fitOverride', 'var') && isfield(fitOverride, 'OneRunNumruns')
        fitS.OneRun.numruns = fitOverride.OneRunNumruns;
    end
    fitS.OneRun.MLNReps = 3; % MLN random restarts for IR fits

    % Set up the replication settings for the chosen weighting
    if fitS.weighting == "depth"
        fitS.OneRun.weightRepDP       = 3; % Decimal places for weight rounding in IR fits
        fitS.OneRun.weightRepInflator = 1; % Weight inflator for IR fits
    elseif fitS.weighting == "age"
        fitS.OneRun.weightRepDP       = 3;
        fitS.OneRun.weightRepInflator = 0.001;
    else
        fitS.OneRun.weightRepDP       = 1;
        fitS.OneRun.weightRepInflator = 1;
    end
end

% --- Chi-squared goodness-of-fit test (optional diagnostic)
fitS.run_chi2gof          = false; % If true, runs chi-squared GOF test on each fit
fitS.chi2binN             = 5;     % Initial number of bins for chi-squared test
fitS.dispChi2             = false; % If true, prints chi-squared results to the console
fitS.chi2MinCountNum      = 5;     % Minimum allowed count per bin (bins are merged if below this)
fitS.chi2MaxCountNum      = 50;    % Maximum allowed count per bin (bins are split if above this)
fitS.enforceBinSizeLimits = true;  % If true, enforces the min/max bin count limits above

% --- Core subset for like-for-like comparisons
fitS.commonCoresMethods = strings(1, 0); % Methods that must all have results for a core to be
                                         %   fitted, e.g. "rejection"; a core any of them skipped
                                         %   is dropped from every fit in this file. Empty = no
                                         %   restriction
fitS.coreFraction       = 1;             % Fraction of the chosen cores to keep, drawn at random
                                         %   (1 = all)
fitS.coreFractionSeed   = 5;             % rng seed for that draw
fitS.chosenCoreNames    = strings(1, 0); % If not empty, fit exactly these cores (after the filters
                                         %   above) - used to give a second file (e.g. the
                                         %   rejection-only noise run) the same cores

% --- Overrides from a driver script
% Overrides of the method lists and core subset from a driver script (e.g. a
% rejection-only run that fits just the rejection methods). No effect
% standalone
if exist('fitOverride', 'var')
    if isfield(fitOverride, 'commonCoresMethods')
        fitS.commonCoresMethods = fitOverride.commonCoresMethods;
    end
    if isfield(fitOverride, 'coreFraction')
        fitS.coreFraction = fitOverride.coreFraction;
    end
    if isfield(fitOverride, 'coreFractionSeed')
        fitS.coreFractionSeed = fitOverride.coreFractionSeed;
    end
    if isfield(fitOverride, 'chosenCoreNames')
        fitS.chosenCoreNames = fitOverride.chosenCoreNames;
    end
    if isfield(fitOverride, 'lin2014Methods')
        fitS.lin2014Methods = fitOverride.lin2014Methods;
    end
    if isfield(fitOverride, 'unthinnedMethods')
        fitS.unthinnedMethods = fitOverride.unthinnedMethods;
    end
end

% Record the fit settings with the results
d.S1.fitS = fitS;

% Choose whether IRfitdists runs in parallel. This is set after d.S1.fitS is
% recorded, so it is not saved with the fit settings
fitS.useParallelComp = false; % If true, uses parfor in IRfitdists (requires Parallel Computing
                              %   Toolbox)

%% Load BIGMACS reference distribution
% Reads the BIGMACS (Lin2014-based) nSR log-normal and transition matrix for
% comparison with the new fits. The MLN fitted here to BM.hist is used later
% as a reference when plotting

% Read the BIGMACS nSR log-normal
BMlognorm = readtable("BIGMACSdata/lognormal.txt"); % (use x values currently used in BIGMACS)
BM.nSR.x  = BMlognorm.Var1';
BM.nSR.px = BMlognorm.Var2;

% Read the Lin2014 sedimentation rate ratios (BM.hist) and the transition matrix
BM.hist = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt");
BM.hist = BM.hist(:,4);
BM.TM   = readmatrix("BIGMACSdata/transition_parameter.txt");

% Fit an MLN to BM.hist, and convert the BIGMACS nSR pdf to log(nSR)
[MLN_BIGMACS, ~, BM.gmfit] = fitMixLogNorm(BM.hist, BM.nSR.x, 2, 3, 343);
[BM.lnSR.x, BM.lnSR.px]    = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);

%% Select core subset for fitting
% Choose the cores that pass the mean SR, water depth and latitude filters
meanSRRestrictionLog   = dataT.meanSR > fitS.minMeanSR;
depthRestrictionLog    = dataT.depths > fitS.minDepth;
latitudeRestrictionLog = abs(dataT.lats) < fitS.maxLat;
d.S1.chooseLog         = meanSRRestrictionLog & depthRestrictionLog & latitudeRestrictionLog;

% Keep only cores that every method in fitS.commonCoresMethods produced
% results for, so that comparisons between methods are on the same cores
for iMethod = 1:numel(fitS.commonCoresMethods)
    % Find the chosen cores this method has no results for
    methodName    = fitS.commonCoresMethods(iMethod);
    hasResultsLog = ~cellfun(@isempty, dataT.(methodName));
    droppedNames  = string(dataT.cores(d.S1.chooseLog & ~hasResultsLog));
    if ~isempty(droppedNames)
        disp("Dropping " + num2str(numel(droppedNames)) + " core(s) with no " + methodName + ...
            " results: " + strjoin(droppedNames, ", "))
    end

    % Remove the cores with no results for this method
    d.S1.chooseLog = d.S1.chooseLog & hasResultsLog;
end

% Keep a random fraction of the chosen cores if asked
if fitS.coreFraction < 1
    % Draw the cores to keep, with fitS.coreFractionSeed so the draw repeats
    chosenIdx = find(d.S1.chooseLog);
    numToKeep = round(fitS.coreFraction.*numel(chosenIdx));
    rng(fitS.coreFractionSeed)
    keptIdx   = sort(randsample(chosenIdx, numToKeep));

    % Keep only those
    d.S1.chooseLog(:)       = false;
    d.S1.chooseLog(keptIdx) = true;
    disp("Keeping a random " + num2str(numToKeep) + " of " + num2str(numel(chosenIdx)) + ...
        " chosen cores")
end

% Fit exactly a given list of cores if asked (so a second file uses the
% same cores as the first)
if ~isempty(fitS.chosenCoreNames)
    d.S1.chooseLog = d.S1.chooseLog & ismember(string(dataT.cores), fitS.chosenCoreNames);
end

% Record which cores were fitted, and how many
d.S1.chosenCoreNames = string(dataT.cores(d.S1.chooseLog))';
d.S1.numCores        = sum(d.S1.chooseLog);
disp("Fitting " + num2str(d.S1.numCores) + " cores")

%% Fit distributions to all nSR estimation methods
% ARfitdists pools all runs across all cores in the chosen subset, applies
% weighting (fitS.weighting), and fits MLN, LN, Gam, and invGam. The Lin2014
% age filter is applied only to the <method>Lin fits; the RSR fits have
% their age differences bounded by the minimum age difference applied while
% sampling, and the rest are fitted as drawn

% Set up the nSR grid on which all fitted PDFs are evaluated
d.logx = -6:0.01:6;
d.x    = exp(d.logx);

% Lin2014 method
% d.S1.LinNSR = ARfitdists(dataT.LinNSR, d.x, d.S1.chooseLog, fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, 1, true, fitS);

%% Set up the number of runs per core for each method (countDivisor)
% ARfitdists divides the number of nSR values by this to get the number of
% nSR calculations per run, used by the fitters as the effective number of
% observations, and weights each core's pairs by countDivisor over that
% core's own number of runs, so cores with different numbers of runs still
% count equally
countDivisors.BMedian   = 1;         % One age model per core
countDivisors.BSamp     = 1000;      % Bchron keeps 1000 iterations per core
countDivisors.rejection = S.numruns;

%% Fit the methods that use the Lin2014 age filter (<method>Lin)
for iMethod = 1:numel(fitS.lin2014Methods)
    % Fit this method with the Lin2014 age filter (applyLin2014AgeFilter = true)
    methodName = fitS.lin2014Methods(iMethod);
    disp("Fitting " + methodName + "Lin")
    d.S1.(methodName + "Lin") = ARfitdists(dataT.(methodName), d.x, d.S1.chooseLog, ...
        fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, ...
        countDivisors.(methodName), true, fitS);
end

%% Fit the sequential RSR methods (RSR<dt><order>)
% No Lin2014 age filter: the minimum age difference applied while sampling
% serves the same purpose. countDivisor is the number of RSR runs per core

% Find every RSR column calcData made
varNames = string(dataT.Properties.VariableNames);
rsrNames = varNames(startsWith(varNames, "RSR"));

% Fit each one
for iMethod = 1:numel(rsrNames)
    methodName = rsrNames(iMethod);
    disp("Fitting " + methodName)
    d.S1.(methodName) = ARfitdists(dataT.(methodName), d.x, d.S1.chooseLog, ...
        fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, S.numruns, false, fitS);
end

%% Fit the methods that sample every age first, as drawn (<method>)
% Every dated depth is used, as in RSR0
for iMethod = 1:numel(fitS.unthinnedMethods)
    methodName = fitS.unthinnedMethods(iMethod);
    disp("Fitting " + methodName)
    d.S1.(methodName) = ARfitdists(dataT.(methodName), d.x, d.S1.chooseLog, ...
        fitS.PooledRuns.weightRepDP, fitS.PooledRuns.weightRepInflator, ...
        countDivisors.(methodName), false, fitS);
end

%%
if fitS.fitIRs
    %% Individual-run (IR) fits - captures run-to-run variability
    % IRfitdists randomly selects one run per core per iteration and fits
    % distributions to that single-run dataset. Repeating this
    % fitS.OneRun.numruns times gives a distribution of fit parameters whose
    % spread quantifies the uncertainty from using any single age model
    % realisation. Done for BSamp (with the Lin2014 age filter) and every RSR
    % column

    % Set up the number of runs, the log(nSR) histogram bins and the methods
    numRuns     = fitS.OneRun.numruns;
    binWidth    = 0.1;             % Bin width in log(nSR) space
    logBinEdges = -5:binWidth:5;   % Log-space bin edges
    irNames     = ["BSamp", rsrNames];

    for iMethod = 1:numel(irNames)
        % Set up the fit name, with the Lin2014 age filter for BSamp only
        methodName = irNames(iMethod);
        if methodName == "BSamp"
            fitName      = "BSampLinIR";
            applyLin2014 = true;
        else
            fitName      = methodName + "IR";
            applyLin2014 = false;
        end

        % Fit each individual run
        disp("Fitting " + fitName)
        rng(2) % Fixed seed for reproducibility
        d.S1.(fitName) = IRfitdists(dataT.(methodName), d.S1.chooseLog, numRuns, d.x, ...
            applyLin2014, fitS);

        % Calculate log-space histogram counts for each individual run
        % (diagnostic), stored as numRuns x numBins; used to visualise the
        % spread of empirical distributions across runs alongside the fitted
        % PDFs
        d.S1.(fitName).lnSRHistCounts = NaN(numRuns, length(logBinEdges)-1);
        for i = 1:numRuns
            run_weightedC  = d.S1.(fitName).weightedC{i};
            run_histCounts = histcounts(log(run_weightedC), 'BinEdges', logBinEdges);
            d.S1.(fitName).lnSRHistCounts(i,:) = run_histCounts ./ ...
                d.S1.fitS.OneRun.weightRepInflator;
        end
    end
end

%% Save all fit results
% Save the d struct (containing all fit outputs, dataT, and settings) to
% fileName2save. The file is not overwritten if it already exists - rename
% fileName2save if you want to rerun and save fresh results
if isfile(fileName2save)
    disp("Warning, not saving file because it already exists. Please rename file")
else
    save(fileName2save, "d", "-v7.3")
end
