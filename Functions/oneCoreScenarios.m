function [core_invSRvals, core_invSRprobs, meanSR, MSI_byage, MSI_bydepth, sedLength, ...
    numDatePairs, ageModes, scenarios, label, numReversals, scenario_meanSR] ...
    = oneCoreScenarios(corename, dataLoc, LabIDs, incDepths, excLabIDs, excDepths, S, plotfigs)
% This function builds the age-depth scenarios for one core and works out
% its summary numbers. It reads and filters the core's radiocarbon dates,
% makes one scenario for each way of choosing a single date at each
% doubly-dated depth (DDD), calibrates the dates, and removes dates from or
% splits any scenario with an age reversal until no scenario has one. It
% then averages each scenario's mean SR, mean sampling interval (MSI),
% sediment length and number of date pairs over the scenarios. Raw 14C ages
% and errors are in kyr; the calibrated age grid and the age modes in yr;
% depths and sediment lengths in cm; SR in cm/kyr; MSI_byage in yr per date
% and MSI_bydepth in cm per date.
%
% dataLoc is "WA" (World Atlas) or "Lin2014"; LabIDs, incDepths, excLabIDs
% and excDepths (from extract3) are passed to filtering; plotfigs turns on
% the diagnostic plots in scenariosDealWithReversals. scenarios holds one
% string array of LabIDs per scenario, ageModes the mode age of each date in
% each scenario, and scenario_meanSR the mean SR of each scenario; oneCoreRSR
% uses all three. label is the LabIDs after
% filtering, with fake labels where they were missing or duplicated.
% core_invSRvals and core_invSRprobs are the scenarios' inverse-SR PDFs
% combined into one (empty unless S.pdfMethod), and numReversals is a
% placeholder (999, not yet implemented). If filtering rejects the core,
% every output except label is empty or NaN; if scenariosDealWithReversals
% drops the core for too many reversals, scenarios is empty and the
% averages are NaN.
%
% S fields used: pdfMethod (S is also passed on to getDataWA, getDatatxt,
% filtering, multiMatcalQ and scenariosDealWithReversals)

%% Read in radiocarbon data
% Read in the radiocarbon ages, depths, errors and LabIDs for this core
% (from the World Atlas NetCDF file or the Lin2014 text file)
if dataLoc == "WA"
    [age_kyr, depth_cm, ageErr_kyr, label] = getDataWA(corename, S);
elseif dataLoc == "Lin2014"
    [age_kyr, depth_cm, ageErr_kyr, label] = getDatatxt(corename, S);
end

%% Filter ages
% Keep the dates chosen for their material (LabIDs, incDepths), remove the
% manually excluded ones (excLabIDs, excDepths), only keep dates within
% S.c14AgeLim, and only keep cores with S.minNumberOfAges or more accepted
% dates
[age_kyr, depth_cm, ageErr_kyr, label, noIncludeInfo, tooFewAges, ~, ~, ~] = filtering(age_kyr, ...
    depth_cm, ageErr_kyr, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

% Stop here if filtering of ages leads to core rejection
if noIncludeInfo == 1 || tooFewAges == 1
    core_invSRvals  = [];
    core_invSRprobs = [];
    meanSR          = NaN;
    MSI_byage       = NaN;
    MSI_bydepth     = NaN;
    numDatePairs    = NaN;
    ageModes        = [];
    sedLength       = NaN;
    scenarios       = [];
    numReversals    = NaN;
    scenario_meanSR = [];
    return
end

%% Create fake labels if needed
% If there are no LabIDs, create temporary LabIDs so that scenarios can be
% made (oneCoreRSR copies this block and relies on the same numbering, so
% change both together)
% FIXME: label is a column, so the for loop below runs only once, with all
% the indices at once, and num2str pads them to one width ("FakeLabel 1"
% ... "FakeLabel10"). The labels are still unique and match the copies
if sum(contains(label, "NaN")) ~= 0
    for iLabel = find(contains(label, "NaN"))
        label(iLabel) = "FakeLabel" + num2str(iLabel);
    end
end

% Check to make sure there aren't duplicated labels in the core data. If
% there are, replace every label with a numbered fake one so they don't
% mess up tests for uniqueness later on
[~, uniqueIdx] = unique(label);
if length(uniqueIdx) ~= length(label)
    for iLabel = 1:length(label)
        label(iLabel) = "xFakeLabel" + num2str(iLabel);
    end
end

%% Create scenarios to get around doubly-dated depths
% Make one scenario (a string array of LabIDs) for each way of choosing a
% single date at each DDD
[scenarios, duplicatedDepths_cm, chosenLabels] = scenariosDDD(depth_cm, label, corename);

%% Calibrate the ages
% Calibrate every date (one probability column per date, on the calAge_yr
% grid)
dateIdx = 1:length(age_kyr);
[ageProbAll, calAge_yr] = multiMatcalQ(age_kyr, ageErr_kyr, dateIdx, S);

%% Set up the scenario struct
% Bundle all per-scenario information into a single struct (SS) so it can
% be passed to removeDuplicateScenarios and scenariosDealWithReversals
numScenarios    = length(scenarios);
SS.scenarios    = scenarios;
SS.CFR          = zeros(numScenarios, 1);   % 1 once a scenario is confirmed reversal-free
SS.chosenLabels = chosenLabels;
SS.invSRvals    = cell(numScenarios, 1);
SS.invSRprobs   = cell(numScenarios, 1);
SS.aveSR        = nan(numScenarios, 1);
SS.numdatepairs = nan(numScenarios, 1);
SS.ageModes     = cell(numScenarios, 1);
SS.lengthSed    = nan(numScenarios, 1);
SS.MSI_byage    = nan(numScenarios, 1);
SS.MSI_bydepth  = nan(numScenarios, 1);

% Initialise the caches of date pairs (by LabID) and their age-difference
% PDFs, so that a pair shared by several scenarios is only worked out once
IDpairs  = strings(0);
agediffV = {};

%% Remove or split scenarios with age reversals
% scenariosDealWithReversals checks the scenarios one at a time. At the
% first reversal it finds, it removes a date from that scenario or splits it
% into new scenarios and returns newScenIndicator = 1, so loop until every
% scenario has been confirmed reversal-free (or the core has been dropped
% for having too many reversals, which empties SS)
newScenIndicator = 1;
while newScenIndicator == 1
    % Remove any scenario that is a copy of another
    SS = removeDuplicateScenarios(SS);

    % Check the scenarios for reversals and deal with the first one found
    [SS, newScenIndicator, IDpairs, agediffV] = ...
        scenariosDealWithReversals(SS, depth_cm, ageProbAll, calAge_yr, ...
        label, corename, duplicatedDepths_cm, IDpairs, agediffV, S, plotfigs);

    % DEBUG: put a breakpoint on the line below to stop when a core is left
    % with no scenarios
    if sum(size(SS.ageModes)) == 0
        a = 1; %#ok<NASGU> % This line used for debugging
    end
end

%% Unpack the scenario struct
% Get each scenario's results out of SS for the calculations below
scenarios             = SS.scenarios;
scenario_meanSR       = SS.aveSR;
scenario_numDatePairs = SS.numdatepairs;
scenario_sedLength    = SS.lengthSed;
scenario_MSI_byage    = SS.MSI_byage;
scenario_MSI_bydepth  = SS.MSI_bydepth;
ageModes              = SS.ageModes;
scenario_invSRvals    = SS.invSRvals;
scenario_invSRprobs   = SS.invSRprobs;

%% Average the summary numbers over the scenarios
% Average each scenario's value, skipping NaNs
meanSR       = mean(scenario_meanSR, 'omitmissing');
numDatePairs = mean(scenario_numDatePairs, 'omitmissing');
sedLength    = mean(scenario_sedLength, 'omitmissing');
numReversals = 999;   % Placeholder, not yet implemented
MSI_byage    = mean(scenario_MSI_byage, 'omitmissing');
MSI_bydepth  = mean(scenario_MSI_bydepth, 'omitmissing');

%% Combine the scenarios' inverse-SR PDFs into a single core PDF
% Only for the analytical PDF method (S.pdfMethod); combinepdfs weights each
% scenario by its sediment length
% FIXME: with S.pdfMethod true, scenariopdfNorm (called by
% scenariosDealWithReversals above) reads S.weighting, which calcData does
% not set (weighting is fitS.weighting in fitData), so the core stops with
% an error before it gets here
if S.pdfMethod
    if ~isempty(scenarios)
        [core_invSRvals, core_invSRprobs] = combinepdfs(scenario_invSRvals, ...
            scenario_invSRprobs, scenario_sedLength);
    else
        core_invSRvals  = [];
        core_invSRprobs = [];
    end
else
    core_invSRvals  = [];
    core_invSRprobs = [];
end

end
