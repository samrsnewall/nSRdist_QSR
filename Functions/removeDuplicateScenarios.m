function SS = removeDuplicateScenarios(SS)
% This function removes every scenario in SS that is a copy of another (the
% same LabIDs in the same order), keeping the first of each set of copies,
% and removes the copies' entries from the other SS fields too. It compares
% the scenarios one length (number of dates) at a time. If SS has no
% scenarios it is returned unchanged.
%
% SS fields used (each is cut down to the same scenarios): scenarios, CFR,
% chosenLabels, invSRvals, invSRprobs, aveSR, numdatepairs, ageModes,
% lengthSed, MSI_byage, MSI_bydepth

% Stop here if there are no scenarios (scenariosDealWithReversals empties SS
% when a core has too many reversals)
if isempty(SS.scenarios)
    return
end

%% Find the scenarios to keep
% Group the scenarios by their number of dates
scenario_numDates = cellfun(@length, SS.scenarios);
uniqueNumDates    = unique(scenario_numDates);

% Initialise the indices of the scenarios to keep
scenarios2keepIdx = [];

% Go through each group and keep the first of each set of identical
% scenarios
for i = 1:length(uniqueNumDates)
    % Get the scenarios with this number of dates, and where they are in
    % SS.scenarios (they need not be together: scenariosDealWithReversals
    % can shorten a scenario or add split ones at the end)
    sameLengthLog       = scenario_numDates == uniqueNumDates(i);
    sameLengthIdx       = find(sameLengthLog);
    sameLengthScenarios = SS.scenarios(sameLengthLog);

    % Put them side by side in a string array, one column per scenario
    sameLengthArray = strings(uniqueNumDates(i), sum(sameLengthLog));
    for j = 1:sum(sameLengthLog)
        sameLengthArray(:,j) = sameLengthScenarios{j};
    end

    % Find the first of each set of identical scenarios (rows once
    % transposed), and add its index in SS.scenarios to the list of
    % scenarios to keep
    [~, uniqueIdx, ~] = unique(sameLengthArray', 'rows');
    scenarios2keepIdx = [scenarios2keepIdx; sameLengthIdx(uniqueIdx)]; %#ok<AGROW>
end

% Keep the scenarios in their original order
scenarios2keepIdx = sort(scenarios2keepIdx);

%% Remove the copies from every field of SS
SS.scenarios    = SS.scenarios(scenarios2keepIdx);
SS.CFR          = SS.CFR(scenarios2keepIdx);
SS.chosenLabels = SS.chosenLabels(scenarios2keepIdx);
SS.invSRvals    = SS.invSRvals(scenarios2keepIdx);
SS.invSRprobs   = SS.invSRprobs(scenarios2keepIdx);
SS.aveSR        = SS.aveSR(scenarios2keepIdx);
SS.numdatepairs = SS.numdatepairs(scenarios2keepIdx);
SS.ageModes     = SS.ageModes(scenarios2keepIdx);
SS.lengthSed    = SS.lengthSed(scenarios2keepIdx);
SS.MSI_byage    = SS.MSI_byage(scenarios2keepIdx);
SS.MSI_bydepth  = SS.MSI_bydepth(scenarios2keepIdx);

end
