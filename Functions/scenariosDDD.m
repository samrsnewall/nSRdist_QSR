function [scenarios, duplicatedDepths_cm, chosenLabels] = scenariosDDD(depth_cm, label, corename)
% This function makes one scenario for each way of choosing a single date at
% each doubly-dated depth (DDD) of a core, using scenariomaker. A DDD is
% where a date's depth is 0.00001 cm or less above the next date's (a
% tolerance for floating-point imprecision in the depths), so the dates must
% be in order of ascending depth. Depths in cm.
%
% scenarios holds one string array of LabIDs per scenario, and chosenLabels
% the LabIDs chosen at the DDDs in each scenario (scenariosDealWithReversals
% uses them to tell which date of a reversal pair to remove). If the core
% has no DDDs, there is one scenario holding every date, chosenLabels is
% {[]} and duplicatedDepths_cm is empty. corename is not used.

%% Find the doubly-dated depths
% Find each pair of neighbouring dates whose depths are the same (to within
% 0.00001 cm), then get the depths that have more than one date
% FIXME: ismember below needs an exact match, so two dates less than
% 0.00001 cm apart but not at exactly the same depth are not grouped (both
% stay in every scenario), and a negative dx_cm (depths not in ascending
% order) is also counted as a DDD
dx_cm               = diff(depth_cm);
sameDepthPairLog    = dx_cm <= 0.00001;           % One entry per pair of neighbouring dates
sameDepthDateLog    = [sameDepthPairLog' false];  % One entry per date (the last has no next date)
duplicatedDepths_cm = unique(depth_cm(sameDepthDateLog));

%% Make the scenarios
if ~isempty(duplicatedDepths_cm)
    % Collect the LabIDs of the dates at each DDD
    DDDlabels = cell(1, length(duplicatedDepths_cm));
    for i = 1:length(duplicatedDepths_cm)
        atThisDepthLog = ismember(depth_cm, duplicatedDepths_cm(i));
        DDDlabels{i}   = label(atThisDepthLog);
    end

    % Make one scenario for each way of choosing one date at each DDD
    [scenarios, chosenLabels, ~] = scenariomaker(DDDlabels, [], label);

else
    % If there are no DDDs, make one scenario holding every date
    scenarios{1}    = label;
    chosenLabels{1} = [];
end

end
