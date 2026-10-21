function[age_kyr, depth_cm, ageErr_kyr, label, noIncludeInfo, tooFewAges, EM, ManE, NotCCR, ageGapTooHigh] = filtering(age_kyr, depth_cm, ageErr_kyr, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S)
% This function filters the radiocarbon dates of one core. It keeps the
% dates chosen for their material (the LabIDs or incDepths listed for the
% core, or every date if LabIDs is "all"), removes the dates manually
% excluded by LabID or depth (excLabIDs, excDepths) and those whose age plus
% or minus its error falls outside S.c14AgeLim, and flags the core if fewer
% than S.minNumberOfAges dates remain. Ages and errors are uncalibrated 14C
% kyr BP (as getDataWA and getDatatxt return them); depths in cm; incDepths
% and excDepths are in m.
%
% noIncludeInfo is true if the core has neither LabIDs nor incDepths to
% choose dates by, and tooFewAges is true if fewer than S.minNumberOfAges
% dates remain; the callers skip the core if either is true. EM, ManE and
% NotCCR hold the dates (age, depth, error, label) left out for their
% material, by manual exclusion and for being outside S.c14AgeLim, each
% taken from the full set of dates, so a date can be in more than one.
% ageGapTooHigh holds the kept dates either side of a gap of more than
% 5 kyr (14C) between neighbouring ages (noted only, not removed).
%
% S fields used: c14AgeLim, minNumberOfAges

%% Initialise the rejection flags
% Set both flags to false (true means the caller should skip this core)
noIncludeInfo = false;   % Signals that the core has no dates chosen (true = no data)
tooFewAges    = false;   % Signals that the core has fewer than S.minNumberOfAges dates

%% Choose dates desired (i.e. planktonic foraminifera dates only)
% Filter the data so it's only those with MSPF
if strcmp(string(LabIDs), "all")
    % If LabIDs says "all", no dates are excluded for their material
    chosenMaterialLog = true(length(age_kyr), 1);

    % Fill the excluded material (EM) outputs (empty here)
    EM.age   = age_kyr(~chosenMaterialLog);
    EM.depth = depth_cm(~chosenMaterialLog);
    EM.error = ageErr_kyr(~chosenMaterialLog);
    EM.label = label(~chosenMaterialLog);
elseif ~isempty(LabIDs) & ~isnan(LabIDs)
    % If LabIDs is not empty, keep the dates with those LabIDs

    % Get the LabIDs of the dates to keep
    LabIDsSplit_Raw = split(LabIDs, ', ');

    % Remove problematic whitespace that isn't removed by strip (char 65279
    % causes many problems)
    problemWhitespaceLog = contains(LabIDsSplit_Raw, char(65279));
    if sum(problemWhitespaceLog) > 0
        disp("Removing problematic whitespace")
        LabIDsSplit_noProblemWhitespace = cell(length(LabIDsSplit_Raw),1);
        for i = 1:length(LabIDsSplit_Raw)
            LabIDcharProb  = LabIDsSplit_Raw{i};
            problemCharLog = ismember(LabIDcharProb, char(65279));
            LabIDcharFix   = LabIDcharProb(~problemCharLog);
            LabIDsSplit_noProblemWhitespace{i} = LabIDcharFix;
        end
        LabIDsSplitGood = LabIDsSplit_noProblemWhitespace;
    else
        LabIDsSplitGood = LabIDsSplit_Raw;
    end

    % Find the dates whose labels are in the LabIDs field
    labelStripped     = strip(string(label));
    LabIDsStripped    = strip(string(LabIDsSplitGood));
    chosenMaterialLog = ismember(labelStripped, LabIDsStripped);
    if sum(chosenMaterialLog) ~= length(LabIDsSplitGood)
        % FIXME: the message lacks a space between "core" and the core name
        warning("Some of the LabIDs listed to be chosen in the COPYcorechoices_MSPF file for core" ...
            + corename + " do not match with the LabIDs read in from the WA2022")
    end

    % Fill the excluded material (EM) outputs
    EM.age   = age_kyr(~chosenMaterialLog);
    EM.depth = depth_cm(~chosenMaterialLog);
    EM.error = ageErr_kyr(~chosenMaterialLog);
    EM.label = label(~chosenMaterialLog);
elseif ~isempty(incDepths) & ~isnan(incDepths)
    % If there are depths to include, keep the dates at those depths

    % Choose what depths to include (incDepths in m, compared in cm)
    incDepthsStr = string(split(incDepths, ', '));
    incDepths_m  = double(incDepthsStr);

    % Find the dates at those depths (tolerance of 1e-6 due to low precision
    % in core GEOFARK12)
    chosenMaterialLog = ismembertol(depth_cm, incDepths_m.*100, 1e-6);

    % Fill the excluded material (EM) outputs
    EM.age   = age_kyr(~chosenMaterialLog);
    EM.depth = depth_cm(~chosenMaterialLog);
    EM.error = ageErr_kyr(~chosenMaterialLog);
    EM.label = label(~chosenMaterialLog);
else
    % If there is no information to include, signal the caller to break out
    % of its function (no dates are chosen, so all go into EM)
    noIncludeInfo     = true;
    chosenMaterialLog = false(size(age_kyr));
    EM.age   = age_kyr;
    EM.depth = depth_cm;
    EM.error = ageErr_kyr;
    EM.label = label;
end

%% Remove manually excluded (ManE) dates
% Initialise the logicals for dates excluded by LabID and by depth
excludedByLabIDLog = false(length(age_kyr), 1);
excludedByDepthLog = false(length(age_kyr), 1);

% Remove the LabIDs of ages that are clearly erroneous (large age
% reversal - manually chosen)

% --- Specified by LabID
if ~isempty(excLabIDs) & ~isnan(excLabIDs)
    % FIXME: contains matches substrings, so an excluded LabID that is part
    % of a longer LabID (e.g. "OS-123" in "OS-1234") excludes that date too
    excLabIDsSplit     = string(split(excLabIDs, ', '));
    excludedByLabIDLog = contains(string(label), excLabIDsSplit);
end

% --- Specified by depth (excDepths in m, compared in cm)
if ~isempty(excDepths) & ~isnan(excDepths)
    if isa(excDepths, 'char')
        excDepthsStr       = string(split(excDepths, ', '));
        excDepths_m        = double(excDepthsStr);
        excludedByDepthLog = ismembertol(depth_cm, excDepths_m.*100, 0.00001);
    elseif isa(excDepths, 'double')
        excDepths_m        = excDepths;
        excludedByDepthLog = ismembertol(depth_cm, excDepths_m.*100, 0.00001);
    else
        % Any other type: carry on with no dates excluded by depth
        disp("excDepths data is not being used for a core because of data type - see filtering function")
    end
end

% Combine the two logicals and set these into Manually Excluded category
manuallyExcludedLog = excludedByLabIDLog == 1 | excludedByDepthLog == 1;
ManE.age   = age_kyr(manuallyExcludedLog);
ManE.depth = depth_cm(manuallyExcludedLog);
ManE.error = ageErr_kyr(manuallyExcludedLog);
ManE.label = label(manuallyExcludedLog);

%% Filter ages to fit within calibration curve
% Keep only the dates whose age minus its error and age plus its error are
% both within S.c14AgeLim (14C kyr BP), to ensure they fit within the limits
% of the calibration curve
youngestAge_kyr = age_kyr-ageErr_kyr;
oldestAge_kyr   = age_kyr+ageErr_kyr;
inCalCurveLog   = youngestAge_kyr >= S.c14AgeLim(1) & oldestAge_kyr <= S.c14AgeLim(2);

% Fill the outside calibration curve range (NotCCR) outputs
NotCCR.depth = depth_cm(~inCalCurveLog);
NotCCR.age   = age_kyr(~inCalCurveLog);
NotCCR.error = ageErr_kyr(~inCalCurveLog);
NotCCR.label = label(~inCalCurveLog);

%% Choose dates that fit all criteria
% Keep the dates chosen for their material, not manually excluded and within
% S.c14AgeLim
chosenLog = chosenMaterialLog & ~manuallyExcludedLog & inCalCurveLog;

% Fill outputs
age_kyr    = age_kyr(chosenLog);
depth_cm   = depth_cm(chosenLog);
ageErr_kyr = ageErr_kyr(chosenLog);
label      = label(chosenLog);

%% Note where there are age differences greater than 5 kyr
% Uses the mean uncalibrated radiocarbon ages, so the 5 kyr threshold is in
% 14C kyr; these dates are noted in ageGapTooHigh, not removed
% TODO: move the 5 kyr age-gap threshold into S

% Calculate differences between means of uncalibrated radiocarbon age (kyr)
dt_kyr = diff(age_kyr);

% Find the pairs of neighbouring ages more than 5 kyr apart
largeGapPairLog = dt_kyr > 5;

% Initiate a logical (default = false) to denote which ages are involved in
% the differences greater than 5 kyr, and mark both ages of each such pair
largeGapAgeLog = false(size(age_kyr));
if sum(largeGapPairLog) ~= 0
    for i = 1:length(dt_kyr)
        if largeGapPairLog(i) == 1
            largeGapAgeLog(i)   = true;
            largeGapAgeLog(i+1) = true;
        end
    end
end

% Fill the large age gap (ageGapTooHigh) outputs
ageGapTooHigh.depth = depth_cm(largeGapAgeLog);
ageGapTooHigh.age   = age_kyr(largeGapAgeLog);
ageGapTooHigh.error = ageErr_kyr(largeGapAgeLog);
ageGapTooHigh.label = label(largeGapAgeLog);

%% Break the core at place of large age gap and see if portion can be used?

%% Filter cores to have at least S.minNumberOfAges dates remaining
% (arbitrary choice... chosen to exclude cores with v low data)
% Flag cores that have fewer than S.minNumberOfAges dates remaining
if length(age_kyr) < S.minNumberOfAges
    tooFewAges = true;
end

end
