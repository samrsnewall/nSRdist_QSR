function [SS, newScenIndicator, IDpairs, agediffV] = scenariosDealWithReversals(SS, depth_cm, ...
    ageProbAll, calAge_yr, label, corename, duplicatedDepths_cm, IDpairs, agediffV, S, plotfigs)
% This function checks the scenarios in SS for age reversals, one at a time
% with scenariopdfNorm, and stores each reversal-free scenario's results in
% SS with SS.CFR = 1 (checked for reversals). At the first scenario with a
% reversal it removes a date from that scenario (when one date of a
% reversal pair is at a doubly-dated depth, DDD) or splits it into new
% scenarios with scenariomaker, and returns straight away with
% newScenIndicator = 1, so that oneCoreScenarios calls it again. If that
% scenario has more than 12 reversals it empties SS instead (the core is
% dropped) and returns newScenIndicator = 0, as it does once every scenario
% has SS.CFR = 1.
%
% Depths and sediment lengths (lengthSed) in cm; calibrated ages and age
% modes in yr; aveSR in cm/kyr; MSI_byage in yr per date and MSI_bydepth in
% cm per date. invSRvals and invSRprobs are empty unless S.pdfMethod.
% IDpairs and agediffV are the caches of date pairs (by LabID) and their
% age-difference PDFs, passed on to scenariopdfNorm (and agediffcalc) and
% returned with any new pairs added; plotfigs turns on scenariopdfNorm's
% diagnostic plots.
%
% SS fields used and set: scenarios, CFR, chosenLabels, invSRvals,
% invSRprobs, aveSR, numdatepairs, ageModes, lengthSed, MSI_byage,
% MSI_bydepth
% S fields used: none directly (S is passed on to scenariopdfNorm, which
% reads reversalCriteria, pdfMethod and weighting, and pdfMinVal through
% agediffcalc)

%% Unpack the scenario struct
% Get each scenario's information out of SS (one entry per scenario)
scenarios    = SS.scenarios;
scenariosCFR = SS.CFR;
chosenLabels = SS.chosenLabels;
invSRvals    = SS.invSRvals;
invSRprobs   = SS.invSRprobs;
aveSR        = SS.aveSR;
numDatePairs = SS.numdatepairs;
ageModes     = SS.ageModes;
sedLength    = SS.lengthSed;
MSI_byage    = SS.MSI_byage;
MSI_bydepth  = SS.MSI_bydepth;

%% Initialise the checked scenario information
% One entry per scenario, filled in as each scenario is checked (or copied
% across if it was confirmed reversal-free in an earlier call)
numScenarios    = length(scenarios);
invSRvalsNew    = cell(numScenarios, 1);
invSRprobsNew   = cell(numScenarios, 1);
scenario_labels = cell(numScenarios, 1);
aveSRNew        = nan(numScenarios, 1);
numDatePairsNew = nan(numScenarios, 1);
ageModesNew     = cell(numScenarios, 1);
sedLengthNew    = nan(numScenarios, 1);
ageLengthNew    = nan(numScenarios, 1); %#ok<NASGU> % Age span of each scenario (yr), not used
MSI_bydepthNew  = nan(numScenarios, 1);
MSI_byageNew    = nan(numScenarios, 1);
scenariosNew    = cell(numScenarios, 1);
scenariosNewCFR = zeros(numScenarios, 1);
chosenLabelsNew = cell(numScenarios, 1);

%% Check each scenario for reversals
% Go through the scenarios in order, and deal with the first reversal found
for i = 1:numScenarios
    if scenariosCFR(i) == 0
        % --- Check a scenario not yet confirmed reversal-free
        % Find this scenario's dates, and get their LabIDs in the order of
        % label
        inScenarioLog      = ismember(string(label), scenarios{i});
        scenarioDateIdx    = find(inScenarioLog == 1);
        scenario_labels{i} = label(inScenarioLog == 1);

        % Choose ages in this scenario
        ageProb = ageProbAll(:,scenarioDateIdx);

        % Check to make sure there aren't duplicated labels in the scenario. If
        % there are, remove one (oneCoreScenarios makes every label unique, so
        % this does not happen)
        [uniqueLabels, uniqueIdx] = unique(scenario_labels{i});
        if length(uniqueIdx) ~= length(scenario_labels{i})
            scenario_labels{i} = uniqueLabels;
        end

        % Calculate the scenario's inverse-SR PDF and summary numbers, and
        % find which pairs of neighbouring dates are reversals
        [invSRvalsNew{i}, invSRprobsNew{i}, aveSRNew(i), reversalPairs, ...
            numDatePairsNew(i), ageModesNew{i}, sedLengthNew(i), ageLengthNew(i), ...
            MSI_byageNew(i), MSI_bydepthNew(i), IDpairs, agediffV] = ...
            scenariopdfNorm(depth_cm, scenarioDateIdx, label, ageProb, calAge_yr, ...
            IDpairs, agediffV, S, plotfigs);

        % DEBUG: put a breakpoint on the line below to stop when a scenario
        % has no age modes
        if sum(size(ageModesNew{i})) == 0
            a = 1; %#ok<NASGU>
        end

        % Find number of reversals
        numReversals = sum(reversalPairs);

        % --- Drop the core if this scenario has more than 12 reversals
        % Empty every field of SS and return newScenIndicator = 0, which ends
        % the loop in oneCoreScenarios
        % FIXME: the message says scenarios where it means reversals, and has
        % no space before the core name
        if numReversals > 12
            disp("Too many scenarios in core" + corename + ": " + num2str(numReversals))
            SS.scenarios    = cell([]);
            SS.CFR          = [];
            SS.chosenLabels = cell([]);
            SS.invSRvals    = cell([]);
            SS.invSRprobs   = cell([]);
            SS.aveSR        = nan(0,1);
            SS.numdatepairs = nan(0,1);
            SS.ageModes     = [];
            SS.lengthSed    = nan(0,1);
            SS.MSI_byage    = nan(0,1);
            SS.MSI_bydepth  = nan(0,1);
            newScenIndicator = 0;
            return
        end

        % --- Remove a date from, or split, a scenario with reversals
        if numReversals ~= 0
            % Note down the LabIDs of the two dates of each reversal pair (one
            % row per pair)
            reversalLabels = strings(numReversals, 2);
            for j = 1:numReversals
                reversalIdx         = find(reversalPairs == 1);
                reversalLabels(j,1) = scenario_labels{i}(reversalIdx(j));
                reversalLabels(j,2) = scenario_labels{i}(reversalIdx(j)+1);
            end

            % Find out if any of the reversal pairs are related to
            % doubly-dated depths (DDDs), using the LabIDs of the dates at
            % depths with only one date
            uniqueDepthLabels = label(~ismember(depth_cm, duplicatedDepths_cm));
            if ~isempty(duplicatedDepths_cm)
                % Find the reversal pairs with one date at a DDD and the other
                % at a depth with only one date
                atUniqueDepthLog = ismember(string(reversalLabels), uniqueDepthLabels);
                oneDateAtDDDLog  = sum(atUniqueDepthLog, 2) == 1;

                if ismember(1, oneDateAtDDDLog)
                    % One date in the reversal pair is from a doubly-dated
                    % depth: remove the other date from the scenario (the one
                    % that is not among the scenario's chosenLabels, the dates
                    % chosen at the DDDs)
                    reversalLabels = reversalLabels(oneDateAtDDDLog, :);
                    date2removeLog = ~ismember(reversalLabels, chosenLabels{i});
                    date2remove    = reversalLabels(date2removeLog);
                    dates2keepLog  = ~ismember(scenarios{i}, date2remove);
                    scenarios{i}   = scenarios{i}(dates2keepLog);

                    % Put the shortened scenario back in SS to be checked
                    % again (CFR 0), keeping the results of the scenarios
                    % before it, and return
                    SS.scenarios    = scenarios;
                    SS.CFR          = [scenariosNewCFR(1:i-1); scenariosCFR(i:end)];
                    SS.chosenLabels = chosenLabels;
                    SS.invSRvals    = invSRvalsNew;
                    SS.invSRprobs   = invSRprobsNew;
                    SS.aveSR        = aveSRNew;
                    SS.numdatepairs = numDatePairsNew;
                    SS.ageModes     = ageModesNew;
                    SS.lengthSed    = sedLengthNew;
                    SS.MSI_byage    = MSI_byageNew;
                    SS.MSI_bydepth  = MSI_bydepthNew;
                    newScenIndicator = 1;
                    return
                end
                % If no reversal pair has exactly one date at a DDD (both dates
                % are at DDDs, or neither is), split the scenario below
            end

            % --- Split the scenario at its reversal pairs
            % Make one new scenario for each way of keeping one date of each
            % reversal pair (scenariomaker)
            reversalPairLabels = cell(1, numReversals);
            for j = 1:sum(reversalPairs)
                reversalPairLabels{1,j} = reversalLabels(j,:)';
            end
            [splitScenarios, ~, ~] = scenariomaker([], reversalPairLabels, scenario_labels{i});
            numSplitScenarios      = length(splitScenarios);

            % Give each new scenario this scenario's chosenLabels
            splitChosenLabels(1:numSplitScenarios,1) = chosenLabels(i);

            % Put the scenarios before this one (checked in this call or an
            % earlier one, so everything about them is in the ...New arrays),
            % the ones after it (not checked yet in this call, so still as
            % they came in) and the new ones (to be checked, CFR 0) together
            % in SS, and return
            SS.scenarios    = [scenariosNew(1:i-1);    scenarios(i+1:end);    splitScenarios];
            SS.CFR          = [scenariosNewCFR(1:i-1); scenariosCFR(i+1:end); zeros(size(splitScenarios))];
            SS.chosenLabels = [chosenLabelsNew(1:i-1); chosenLabels(i+1:end); splitChosenLabels];
            SS.invSRvals    = [invSRvalsNew(1:i-1);    invSRvals(i+1:end);    cell(numSplitScenarios,1)];
            SS.invSRprobs   = [invSRprobsNew(1:i-1);   invSRprobs(i+1:end);   cell(numSplitScenarios,1)];
            SS.aveSR        = [aveSRNew(1:i-1);        aveSR(i+1:end);        nan(numSplitScenarios,1)];
            SS.numdatepairs = [numDatePairsNew(1:i-1); numDatePairs(i+1:end); nan(numSplitScenarios,1)];
            SS.ageModes     = [ageModesNew(1:i-1);     ageModes(i+1:end);     cell(numSplitScenarios,1)];
            SS.lengthSed    = [sedLengthNew(1:i-1);    sedLength(i+1:end);    nan(numSplitScenarios,1)];
            SS.MSI_byage    = [MSI_byageNew(1:i-1);    MSI_byage(i+1:end);    nan(numSplitScenarios,1)];
            SS.MSI_bydepth  = [MSI_bydepthNew(1:i-1);  MSI_bydepth(i+1:end);  nan(numSplitScenarios,1)];
            newScenIndicator = 1;
            return
        end

        % --- Confirm a scenario with no reversals
        % Its results are already in the ...New arrays, so mark it as
        % confirmed reversal-free
        scenariosNew(i)    = scenarios(i);
        scenariosNewCFR(i) = 1;
        chosenLabelsNew(i) = chosenLabels(i);

    else
        % --- Copy a scenario confirmed reversal-free in an earlier call
        % Copy all its information straight through
        scenariosNew(i)    = scenarios(i);
        scenariosNewCFR(i) = 1;
        ageModesNew(i)     = ageModes(i);
        chosenLabelsNew(i) = chosenLabels(i);
        invSRvalsNew{i}    = invSRvals{i};
        invSRprobsNew{i}   = invSRprobs{i};
        aveSRNew(i)        = aveSR(i);
        numDatePairsNew(i) = numDatePairs(i);
        sedLengthNew(i)    = sedLength(i);
        MSI_bydepthNew(i)  = MSI_bydepth(i);
        MSI_byageNew(i)    = MSI_byage(i);
    end

    % No reversal found so far
    newScenIndicator = 0;
end

%% Put the checked scenarios back in SS
% Only reached when every scenario is confirmed reversal-free
SS.scenarios    = scenariosNew;
SS.CFR          = scenariosNewCFR;
SS.chosenLabels = chosenLabelsNew;
SS.invSRvals    = invSRvalsNew;
SS.invSRprobs   = invSRprobsNew;
SS.aveSR        = aveSRNew;
SS.numdatepairs = numDatePairsNew;
SS.ageModes     = ageModesNew;
SS.lengthSed    = sedLengthNew;
SS.MSI_byage    = MSI_byageNew;
SS.MSI_bydepth  = MSI_bydepthNew;

end
