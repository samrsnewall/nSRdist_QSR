%% calcData - stage 1 of the nSRdist pipeline
% This script calculates normalised sedimentation rate (nSR) histories for a
% chosen set of marine sediment cores, using several independent methods,
% and saves the results to a .mat file in Results/ for fitData to use. So
% calcData saves dataT, S and rawdataManual to Results/, and fitData then
% loads that file, fits distributions and saves the fit results
%
% The nSR methods are:
%   LinNSR    - Attempts an exact replication of the Lin2014 method
%   BMedian   - Uses the median of the Bchron posterior age at each dated
%               depth to calculate pairwise SRs; the mean SR for
%               normalisation is taken from the full Bchron model. Similar
%               in concept to the Lin et al. (2014) approach, though not an
%               exact replication
%   BSamp     - Like BMedian but samples from the full Bchron posterior age
%               distribution (bchronProb) across S.numruns iterations, so it
%               propagates age uncertainty more fully than BMedian. Results
%               are stored in the bchronProb cell array
%   RSR0, RSR500, RSR1000
%             - Random Sampling from calibrated Radiocarbon (RSR) age PDFs,
%               with a minimum age difference (dt_min = 0, 500 or 1000 yr)
%               between consecutive dated depths, applied while the ages are
%               sampled. Uses oneCoreRSR; the restriction filters out pairs
%               too close in age to give reliable SR estimates. Run once for
%               each sequential sampling order in S.sampleOrders: the
%               shallow-to-deep order used in the manuscript is stored as
%               RSR0, RSR500 and RSR1000, and the other orders (sensitivity
%               tests in the response to reviewers) are stored with the
%               order appended (RSR500d2s, RSR500random)
%   rejection - Joint rejection sampling from the same calibrated PDFs
%               (oneCoreRSR with sampleOrder "rejection"), with no
%               restriction, so no sampling order is involved. A sensitivity
%               test of the sequential sampling, compared with RSR0 on the
%               cores where rejection sampling is feasible
%
% Doubly-dated depths and reversals are dealt with before the RSR nSRs are
% calculated: oneCoreScenarios is called for each core and builds
% "scenarios" (all valid orderings of radiocarbon dates that avoid age
% reversals) by handling:
%   - Doubly-dated depths (DDDs): multiple dates at the same depth
%   - Age reversals: probabilistic reversals are resolved by constructing
%     alternative scenario pathways through the dated sequence
% The resulting scenarios are shared by all RSRx methods and inform the
% per-scenario mean SR used for normalisation
%
% Data are loaded from an Excel spreadsheet (S.sheet) which records, for
% each core:
%   - Core name, location, water depth, ocean basin
%   - Which radiocarbon dates to include (planktonic foraminifera only)
%   - Which dates to exclude (benthic foram, outliers, large gaps)
%   - Flags linking rows to the Mulitza 2022 World Atlas (useLin/usePF)
% The raw radiocarbon measurement data are read from the Mulitza 2022
% dataset (S.WApath) or from the Lin2014 supplementary database
%
% Saved to Results/S.dataOutputFile.mat:
%   dataT         - Table with one row per core containing cores, lats,
%                   longs, depths, ocean, meanSR, ageModes, sedimentlength,
%                   agecoverage, num14cpairs, MSI_byage, MSI_bydepth,
%                   LabIDs, incDepths, excLabIDs, excDepths, dataLoc,
%                   core_invSRvals, core_invSRprobs, corescenarios,
%                   scenario_meanSR, RSR0, RSR500, RSR1000 (shallow to
%                   deep; the same with the order appended for each other
%                   order in S.sampleOrders), rejection, and BMode, BMedian,
%                   BSamp, LinNSR
%   S             - Settings structure used for this run
%   rawdataManual - Subset of the raw metadata table used
%
% The settings are in S, set in the "Create settings structure" section;
% see the comment on each field there
%
% It needs Functions/ (all helper functions, added to the path at startup),
% MatCal (matcal / matcalq for radiocarbon calibration) and R with Bchron
% (called via Rscript for Bchron age-depth modelling)

%% Add folder of necessary functions to path
addpath('Functions')

%% Create settings structure
% All run parameters are stored in S and saved alongside the results, so
% any output file is self-documenting

% --- Paths
% Find the root of the repository from where this script is saved
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    error('Please run calcData.m as a script file, not pasted into the Command Window.');
end
S.sandboxPath = fileparts(scriptPath);            % Root directory of the repository
S.RscriptPath = "/usr/local/bin/Rscript";         % Path to Rscript; find yours with `which Rscript`
                                                  %   in Terminal
S.WApath      = "WA_Foraminiferal_Isotopes_2022"; % Path to the Mulitza 2022 World Atlas dataset
S.sheet       = "DataSheets/datasheet.xlsx";      % Metadata spreadsheet controlling which cores and
                                                  %   dates are used

% --- Output file
stringID         = "QSR";
S.dataOutputFile = "dataT_" + stringID; % Saved to Results/S.dataOutputFile.mat (the manuscript
                                        %   run; results_QSR.m sets this through calcOverride)

% --- Core selection
S.coreSelection = "good"; % "good" = every core passing the depth / latitude / problem-core filters
                          %   below (the production run, which takes hours); or a string array of
                          %   core names to run a subset, e.g. ["A7" "MD95-2043"] runs in seconds
S.maxNumCores   = [];     % If not empty, keep only the first S.maxNumCores of the chosen cores,
                          %   alphabetically by name

% --- Core selection filters
S.problemCores     = "KNR197-3-36GGC";
                           % Cores left out by name, whatever the filters below say; including
                           %   KNR197-3-36GGC often led to runs crashing
S.minimumCoreDepth = 1000; % Minimum water depth (m); excludes shallow-water cores
S.maxAtlanticLatN  = 40;   % Maximum northern latitude for Atlantic cores (degrees)
S.maxOtherLatN     = 40;   % Maximum northern latitude for non-Atlantic cores
S.maxAtlanticLatS  = 40;   % Maximum southern latitude for Atlantic cores (degrees)
S.maxOtherLatS     = 40;   % Maximum southern latitude for non-Atlantic cores

% --- Data source and pre-processing
S.replicateLin2014  = false;  % If true, overrides other settings to match Lin2014 exactly (see
                              %   block below)
S.useLin            = true;   % Include cores from the Lin2014 database
S.modifyLin2014Data = true;   % If true, applies additional reversals / gap removal on Lin2014 dates
S.usePF             = true;   % Include cores added from the Mulitza 2022 dataset (planktonic foram
                              %   dates only)
S.duplicateLinCores = "KF13"; % Lin2014 cores that are the same core as a WA core (KF13 is
                              %   GEOFARKF13, identical dates); removed when both datasets are used,
                              %   unless replicating Lin2014
S.DeltaRError       = 200;    % 1-sigma error (yr) on the marine reservoir age correction (DeltaR);
                              %   propagated through calibration
S.c14AgeLim         = [0 50]; % Accepted radiocarbon age range (kyr); dates outside this window are
                              %   excluded
S.normWithRunAve    = true;   % If true, normalises each RSR run's SRs by that run's own mean SR;
                              %   if false, uses a single shared mean SR across all runs for a core
S.numruns           = 10000;  % Number of Monte Carlo iterations for BSamp and RSRx methods

% --- Bchron settings
S.useBchron              = true;    % If false, skips Bchron approaches entirely
S.BchronFilter           = true;    % If true, removes obvious outliers and large age gaps before
                                    %   running Bchron
S.BchronFolderName       = "All1_RLGtrue_DS0p05_Jun2";
                                    % Subfolder under BchronFolders/ from which to read pre-computed
                                    %   Bchron outputs
S.BchronOutlier          = 0.05;    % OutlierProbs value passed to Bchron (prior probability of each
                                    %   date being an outlier)
S.BchronReversalCriteria = 0.75;    % Fraction of Bchron MCMC samples that must reject an age for it
                                    %   to be flagged as a reversal
S.BchronCalCurve         = "Marine20";
                                    % Calibration curve to use in Bchron ("Marine20" or "Marine09"
                                    %   for Lin2014 replication)
S.BchronDepthSpacing     = 0.05;    % Depth interval (cm) at which Bchron reports age estimates
                                    %   along the core
S.BchronReDo             = false;   % If false, re-uses existing Bchron outputs if available; if
                                    %   true, reruns Bchron for all cores

% --- RSRx method settings
S.minNumberOfAges     = 4;     % Minimum number of accepted dates a core must have to be processed
S.reversalCriteria    = 0.75;  % Fraction of sampled SR estimates between two adjacent dates that
                               %   must be negative to classify the pair as a reversal and trigger
                               %   scenario branching
S.removeLargeGaps     = true;  % If true, manually excludes date pairs with very large age gaps
                               %   before scenario construction
S.pdfMinVal           = 1e-6;  % Probability density cutoff below which calibrated age PDF tails
                               %   are truncated, reducing vector length and accelerating
                               %   calculations
S.sampleOrders        = ["shallow2deep", "deep2shallow", "random"];
                               % Sequential RSR sampling orders to run; each gives its own
                               %   RSR0/RSR500/RSR1000 columns. shallow2deep (the manuscript
                               %   order) is stored without a suffix (RSR500); the others get the
                               %   order appended (RSR500d2s, RSR500random)
S.rsrMinAgeDiffs      = [0 500 1000];
                               % Minimum age differences (yr) the sequential RSR sampler is run
                               %   with, for each order in S.sampleOrders; a value left out skips
                               %   that RSRx block below, so its column is not added to dataT
S.useRejection        = true;  % If true, also runs the joint rejection sampler once, with no
                               %   restriction (rejection)
S.rejectionBatchSize  = 1000;  % Candidate runs the rejection sampler draws per batch
S.rejectionMaxBatches = 20;    % Batches to try before giving up on a core (its calibrated PDFs
                               %   overlap too much for a monotone run) and skipping it
S.pdfMethod           = false; % If true, runs the full analytical PDF method for all date pairs;
                               %   if false, uses PDF only for DDDs and reversals, sampling
                               %   otherwise
S.useModes            = false; % If true, uses the mode of each calibrated age PDF as a point
                               %   estimate rather than sampling (diagnostic / fast-approximate
                               %   mode)

% --- Lin2014 approach
S.Lin2014Method = false; % Whether to calculate nSRs with Lin2014 method

% --- Other settings
S.matcalFast         = true;  % Whether to use my modified version of matcal. Removes need to add
                              %   MatCal to MATLAB (copy stored within repo). Speeds up scenario
                              %   construction for RSRx methods by ~10x
S.plotAgeModes       = false; % Whether to plot all age modes after creation of scenarios
S.plotAgeDepthModels = false; % Whether to plot calibrated ages against depth for each core

% --- Overrides from a driver script
% If a calling script has defined calcOverride (a struct of S fields), those
% fields replace the values above. Absent when run standalone, so no effect
if exist('calcOverride', 'var')
    overrideFields = fieldnames(calcOverride);
    for iField = 1:numel(overrideFields)
        S.(overrideFields{iField}) = calcOverride.(overrideFields{iField});
    end
end

%% Switch to the Lin2014 set-up if replicating Lin2014
if S.replicateLin2014 == 1
    % Set up the data, calibration and Bchron settings as Lin2014 had them
    S.modifyLin2014Data = false;
    S.DeltaRError       = 0;
    S.c14AgeLim         = [0 55];
    S.useBchron         = true;
    S.BchronCalCurve    = "Marine09";
    S.BchronReDo        = false;

    % Set up core choice settings (first two explicitly stated, others
    % implicitly taken from cores they used)
    S.minimumCoreDepth = 1000; % Minimum core depth (mbsf) in m
    S.maxAtlanticLatN  = 40;   % Maximum north latitude of Atlantic cores (explicitly stated
                               %   in paper)
    S.maxOtherLatN     = 52;   % (not explicitly stated in paper, max LatN is 51ish)
    S.maxAtlanticLatS  = 50;   % (not explicitly stated in paper, max LatS is -46ish)
    S.maxOtherLatS     = 50;   % (not explicitly stated in paper, max LatS is -46ish)
end

%% Load metadata from the spreadsheet
% The spreadsheet (S.sheet) has one row per core and records:
%   CoreName, LatitudeDec, LongitudeDec, WaterDepthM, Basin
%   Lin2014     - 1 if the core was used by Lin et al. (2014)
%   Lin2014Keep - 1 if the core should be kept for replication even if not
%                 in WAuse
%   WAuse       - 1 if the core should be used from the Mulitza 2022 dataset
% For each core, LabIDs of included dates (planktonic foram only) and
% excluded dates (benthic foram, outliers, large gaps) are also stored
rawdata = readtable(S.sheet);

% Select relevant information based on what subsets of data have been chosen
if S.useLin && ~S.usePF
    rawdataManual = rawdata(rawdata.Lin2014 == 1, :);
elseif ~S.useLin && S.usePF
    rawdataManual = rawdata(rawdata.WAuse == 1 & isnan(rawdata.Lin2014), :);
elseif S.useLin && S.usePF
    rawdataManual = rawdata(rawdata.WAuse == 1 | rawdata.Lin2014Keep == 1, :);
end

%% Create logical to choose all good cores
% Find the problem cores named in S.problemCores, to leave out
badLog = ismember(string(rawdataManual.CoreName), S.problemCores);

% Remove Lin2014 cores that duplicate a WA core already in the set, unless
% replicating Lin2014 (KF13 is GEOFARKF13; the WA version is kept because its
% spreadsheet row has the materials and MSPF depths)
if S.useLin && S.usePF && ~S.replicateLin2014
    duplicateLinLog = ismember(string(rawdataManual.CoreName), S.duplicateLinCores);
    badLog          = badLog | duplicateLinLog;
end

% Find which cores are in the Atlantic, and the latitude of every core
inAtlanticLog = contains(string(rawdataManual.Basin), "Atlantic");
allLats       = rawdataManual.LatitudeDec;

% Find the cores within the latitude limits for their basin
AtlanticLatLog   = inAtlanticLog & allLats <= S.maxAtlanticLatN & allLats >= -S.maxAtlanticLatS;
OtherBasinLatLog = ~inAtlanticLog & allLats <= S.maxOtherLatN & allLats >= -S.maxOtherLatS;

% Restrict which cores to analyse based on metadata, such as depth or
% latitude
latitudeRestrictionLog = AtlanticLatLog | OtherBasinLatLog;
depthRestrictionLog    = rawdataManual.WaterDepthM >= S.minimumCoreDepth;

% Therefore create logical of all good cores
withinRestrictionsLog = latitudeRestrictionLog & depthRestrictionLog;
goodLog               = badLog == 0 & withinRestrictionsLog == 1;

%% Get material data from Excel
% Choose the cores: every good core, or the ones named in S.coreSelection
if isequal(S.coreSelection, "good")
    chosenCoresLog = goodLog;
else
    chosenCoresLog = ismember(string(rawdataManual.CoreName), S.coreSelection);
end

% Keep only the first S.maxNumCores chosen cores, alphabetically, if asked
if ~isempty(S.maxNumCores)
    chosenNames    = sort(string(rawdataManual.CoreName(chosenCoresLog)));
    chosenNames    = chosenNames(1:min(S.maxNumCores, numel(chosenNames)));
    chosenCoresLog = ismember(string(rawdataManual.CoreName), chosenNames);
end

% Get relevant metadata from Excel Spreadsheet into useful variables
numCores    = sum(chosenCoresLog);
cores       = table2array(rawdataManual(chosenCoresLog, "CoreName")); % Take list of MSPF core names
lats        = table2array(rawdataManual(chosenCoresLog, "LatitudeDec"));
longs       = table2array(rawdataManual(chosenCoresLog, "LongitudeDec"));
depths      = table2array(rawdataManual(chosenCoresLog, "WaterDepthM"));
ocean       = table2array(rawdataManual(chosenCoresLog, "Basin"));

% Get the dates to include and exclude for each core, and where its data are
[LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(rawdataManual, chosenCoresLog, S);

%% Plot calibrated radiocarbon dates against depth
if S.plotAgeDepthModels
    for iPlot = (1:length(cores))
        corePlotCal(cores{iPlot}, LabIDs{iPlot}, incDepths{iPlot}, excLabIDs{iPlot}, ...
            excDepths{iPlot}, dataLoc(iPlot), S)
    end
end

%% Calculate nSRs from the Bchron age-depth models (BMode, BMedian, BSamp)
% nSRBchron runs Bchron for each core (or loads existing results) and
% returns nSR histories in 3-row matrix format for each method:
%   BMode   - nSR from the modal Bchron age at each dated depth
%   BMedian - nSR from the median Bchron age at each dated depth
%   BProb   - nSRs sampled from the Bchron posterior (S.numruns samples)
%   LinNSR  - nSR history matching the Lin2014 mean SR normalisation
% Skipped entirely if S.useBchron is false

% Initialise cells
BMode         = cell(numCores, 1);
BMedian       = cell(numCores, 1);
BSamp         = cell(numCores, 1);
LinNSR        = cell(numCores, 1);
bchronMeanSRs = NaN(numCores, 1);

% Run (or load) Bchron for each core and calculate its nSRs
if S.useBchron
    for i = 1:numCores
        [BMode{i}, BMedian{i}, BSamp{i}, LinNSR{i}, bchronMeanSRs(i)] = nSRBchron(cores{i}, ...
            LabIDs{i}, incDepths{i}, excLabIDs{i}, excDepths{i}, dataLoc(i), S);
    end
end

%% Build scenarios and summary statistics for each core (oneCoreScenarios)
% oneCoreScenarios is called for every core before the RSRx sampling
% begins. It serves two purposes:
%   1. Scenario construction - resolves doubly-dated depths (DDDs) and age
%      reversals by enumerating all valid date orderings ("scenarios"). Each
%      scenario is a specific assignment of ages to depths that produces no
%      reversals. These scenarios are later used by oneCoreRSR
%   2. Summary statistics - calculates the mean SR, mean sampling interval
%      (MSI by age and by depth), total sediment length, number of 14C date
%      pairs, per-depth age modes, and per-scenario mean SRs used for
%      normalisation in the RSRx methods
% It also calculates an inverse-SR PDF (core_invSRvals / core_invSRprobs)
% with the analytical PDF method for diagnostic purposes (stored in dataT)

if S.replicateLin2014
    % Code has problems when ages get too large, because they get to the
    % edge of the vector that is returned by MatCal
    S.c14AgeLim = [0 48];
end

% Initialise variables to hold this information
core_invSRvals  = cell(numCores,1);
core_invSRprobs = cell(numCores,1);
meanSR          = nan(numCores,1);
MSI_byage       = nan(numCores,1);
MSI_bydepth     = nan(numCores,1);
sedimentlength  = nan(numCores,1);
num14cpairs     = nan(numCores,1);
ageModes        = cell(numCores, 1);
transProbs      = nan(3,3,numCores);
corescenarios   = cell(numCores,1);
newLabels       = cell(numCores,1);
numReversals    = nan(numCores,1);
scenario_meanSR = cell(numCores, 1);

% Calculate SR distribution for each core, as well as meanSR and other
% useful information (the last input, 0, turns off oneCoreScenarios' plots)
tic
for i = 1:numCores
    disp("Working on " + cores{i})
    [core_invSRvals{i}, core_invSRprobs{i}, meanSR(i), MSI_byage(i), ...
        MSI_bydepth(i), sedimentlength(i), num14cpairs(i), ageModes{i}, ...
        corescenarios{i}, newLabels{i}, numReversals(i), scenario_meanSR{i}] ...
        = oneCoreScenarios(cores{i}, dataLoc(i), LabIDs{i}, incDepths{i}, excLabIDs{i}, ...
        excDepths{i}, S, 0);
end
toc
disp("Created all scenarios")

%% Create a table that holds all useful metadata
% Calculate how much time each core covers
agecoverage = sedimentlength./meanSR;

% Put the metadata and scenario information in one table, one row per core
dataT = table(cores, lats, longs, depths, ocean, meanSR, ageModes, ...
    sedimentlength, agecoverage, num14cpairs, MSI_byage, MSI_bydepth, ...
    LabIDs, incDepths, excLabIDs, excDepths, dataLoc, core_invSRvals, ...
    core_invSRprobs, corescenarios, scenario_meanSR);

%% Plot figure showing all the ages being used
if S.plotAgeModes
    plotAgeModes(chosenCoresLog, chosenCoresLog, ageModes, cores)
end

%% Calculate nSRs by random sampling from the calibrated radiocarbon PDFs (RSRx)
% oneCoreRSR samples S.numruns times from the calibrated radiocarbon age
% PDFs of each core and calculates pairwise SRs. The output is a 3-row nSR
% matrix (see README "nSR matrix format") holding all sampled nSR values
% and their associated depth/age differences
%
% Three variants are run for each sequential sampling order in
% S.sampleOrders, differing only in the minimum age difference (dt_min)
% required between consecutive dated depths, which is applied while the
% ages are sampled:
%   RSR0    - no restriction (any dt accepted)
%   RSR500  - pairs with dt < 500 yr excluded
%   RSR1000 - pairs with dt < 1000 yr excluded
% Which of the three are run is set by S.rsrMinAgeDiffs.
% The shallow-to-deep order (the manuscript method) is stored as RSR0,
% RSR500 and RSR1000; the other orders are stored with the order appended
% (RSR500d2s, RSR500random). The joint rejection sampler is then run once
% with no restriction (rejection)
%
% The restriction reduces bias from closely-spaced dates whose calibrated
% age PDFs overlap substantially. All three RSRx loops are parfor for speed (so rng seeds
% do not reproduce their draws; the workers use their own random streams)

for iOrder = 1:numel(S.sampleOrders)
    % Get this sampling order
    sampleOrder = S.sampleOrders(iOrder);

    % Set up the short order name used in the dataT column names (none for
    % shallow2deep, the manuscript order)
    orderTag = replace(sampleOrder, ["shallow2deep", "deep2shallow"], ["", "d2s"]);

    % --- RSR0: no minimum age-difference restriction
    if ismember(0, S.rsrMinAgeDiffs)
        % Initialise variables
        disp("Starting RSR0" + orderTag)
        RSR0 = cell(numCores, 1); % Holds RSR0 data

        parfor i = 1:numCores
            disp(cores{i})
            [RSR0{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, ...
                incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, ...
                S, 0, sampleOrder);
        end
        dataT.("RSR0" + orderTag) = RSR0;
    end

    % --- RSR500: minimum age-difference restriction = 500 yr
    if ismember(500, S.rsrMinAgeDiffs)
        % Initialise variables
        disp("Starting RSR500" + orderTag)
        RSR500 = cell(numCores, 1); % Holds RSR500

        parfor i = 1:numCores
            [RSR500{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, ...
                incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, ...
                S, 500, sampleOrder);
        end
        dataT.("RSR500" + orderTag) = RSR500;
    end

    % --- RSR1000: minimum age-difference restriction = 1000 yr
    if ismember(1000, S.rsrMinAgeDiffs)
        % Initialise variables
        disp("Starting RSR1000" + orderTag)
        RSR1000 = cell(numCores, 1); % Holds RSR1000

        parfor i = 1:numCores
            [RSR1000{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, ...
                incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, ...
                S, 1000, sampleOrder);
        end
        dataT.("RSR1000" + orderTag) = RSR1000;
    end
end

% --- rejection: joint rejection sampler, no restriction
% Draws every date at once and keeps only runs whose ages increase with
% depth, so no sampling order is involved. Cores whose calibrated PDFs
% overlap too much for a monotone run to be found within
% S.rejectionMaxBatches batches are skipped (empty cell)
if S.useRejection
    disp("Starting rejection")
    rejection = cell(numCores, 1); % Holds rejection data

    for i = 1:numCores
        [rejection{i}, ~] = oneCoreRSR(cores{i}, dataLoc(i), corescenarios{i}, LabIDs{i}, ...
            incDepths{i}, excLabIDs{i}, excDepths{i}, scenario_meanSR{i}, ageModes{i}, ...
            S, 0, "rejection");
    end
    dataT.rejection = rejection;
end
disp("All RSRx samples finished")

%% Assemble final table and save
% Append the Bchron outputs to dataT (the RSR and rejection columns were
% added above), then save everything (dataT, S, rawdataManual) to
% Results/S.dataOutputFile.mat. This file is the input to fitData
dataT = addvars(dataT, BMode, BMedian, BSamp, LinNSR);

% Save as -v7.3 because dataT is over the 2 GB limit of the default format
% once every sampling order is included
save(fullfile("Results", S.dataOutputFile), "dataT", "S", "rawdataManual", "-v7.3")
