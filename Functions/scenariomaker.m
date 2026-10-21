function [scenarios, chosenLabels, ungroupedLabels] = scenariomaker(DDDlabels, reversalPairLabels, label)
% This function makes one scenario for each way of choosing a single date
% from each group of dates that cannot all be in the same scenario: the
% dates at each doubly-dated depth (DDDlabels, from scenariosDDD) and the
% two dates of each age reversal pair (reversalPairLabels, from
% scenariosDealWithReversals). Either can be []. Each group is a string
% column of LabIDs, and label holds the LabIDs of the whole core (from
% scenariosDDD) or of the scenario being split (from
% scenariosDealWithReversals).
%
% Each scenario is a string array of LabIDs: the dates of label that are in
% no group (ungroupedLabels, the same in every scenario), followed by the
% dates chosen from the groups (chosenLabels, one per group, in group
% order).

%% Put the two kinds of group together
groupLabels = [DDDlabels reversalPairLabels];

%% Find the dates that are in no group
% These are in every scenario, so they only need finding once. A date is
% grouped only if its whole LabID is in a group: one that merely contains a
% grouped LabID (V35-5's Gsac12 contains its DDD LabID Gsac1) is not
allGroupLabels  = vertcat(groupLabels{1,:});
inGroupLog      = ismember(string(label), allGroupLabels);
ungroupedLabels = label(~inGroupLog);

%% Make every combination of one date from each group
% Get the combinations as a table, one row per scenario and one column per
% group
combinationsT    = combinations(groupLabels{:});
combinationsSize = size(combinationsT);
numScenarios     = combinationsSize(1);

% Initialise cells
scenarios    = cell(numScenarios, 1);
chosenLabels = cell(numScenarios, 1);

% Make each scenario from the ungrouped dates and this combination's dates
for i = 1:numScenarios
    % Get this combination's LabIDs out of the table row as a string column,
    % by joining them with "_" and splitting them again
    % FIXME: a LabID that contains "_" would be split into two
    chosenLabelsT      = combinationsT(i,:);
    chosenLabelsRow    = chosenLabelsT{:,:};
    chosenLabelsJoined = join(chosenLabelsRow, '_');
    chosenLabels{i}    = split(chosenLabelsJoined, '_');

    % Put the ungrouped dates first, then the chosen ones
    scenarios{i} = [ungroupedLabels; chosenLabels{i}];
end

end
