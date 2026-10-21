function[transNums, numCSE2x, TM, TMw] = TMcalculationAR(nSRcountsArray, fitS)
% This function calculates a transition matrix (TM) between the three nSR
% states of Lin et al. (2014), contraction (C), steady (S) and expansion
% (E), from a pooled nSR matrix, laid out like BIGMACS's (BM.TM, from
% BIGMACSdata/transition_parameter.txt): row 1 the initial state, rows 2-4
% the probabilities of going from C, S and E (one row each) to C, S and E
% (one column each). Each pair is put in a state by its nSR, the transitions
% between consecutive pairs of the same run are counted, and each count is
% divided by the number of transitions out of its state; TMw does the same
% but counts the steps within a stretch of one state by fitS.weighting.
% nSR is unit-free; dx in cm; dt in yr (converted to kyr for the "age"
% weights).
%
% nSRcountsArray is a 3-row nSR matrix (nSR, dx_cm, dt_yr) with each run
% starting with a NaN column (ARfitdists passes its pooled array from before
% the Lin2014 age filter); the NaN columns get the state 'b', so no
% transition is counted across a run start. transNums is the 3 x 3 number
% of each transition (rows from, columns to, in the order C, S, E), numCSE2x
% the number of transitions out of C, S and E, TM the 4 x 3 transition
% matrix with row 1 left NaN (see the TODO below), and TMw the 4 x 3
% weighted transition matrix with row 1 the fraction of the weighted
% transitions that start in each state. C is nSR >= 1.085 and E is nSR <
% 0.922 throughout. ARfitdists uses only TM.
%
% fitS fields used: weighting

%% Get the nSR values
% Get the nSR of every pair, keeping the NaN run-start columns (they mark
% the breaks between runs)
nSRs = nSRcountsArray(1,:);

%% Calculate TM
% Put each pair into a state by its nSR, with the thresholds from Lin et al.
% (2014). Categories are Steady, Expansion, Contraction (S, E, C)
% TODO: the 0.922 and 1.085 thresholds could move into fitS
stateStr = '';
for i = 1:length(nSRs)
    if nSRs(i) >= 0.922 && nSRs(i) < 1.085
        stateStr(i) = 'S';
    elseif nSRs(i) >= 0 && nSRs(i) < 0.922
        stateStr(i) = 'E';
    elseif nSRs(i) >= 1.085 && nSRs(i) < inf
        stateStr(i) = 'C';
    else
        % If the nSR is not one of those values it'll be a NaN, which marks
        % the start of a run of ages (and the end of the one before)
        stateStr(i) = 'b';
    end
end

% --- Convert this to a transition matrix
% Find the number of times each transition occurs
% TODO: Currently this counts events from shallower to deeper..., needs to
% be the other way around. So flip the whole vector and then do it (needs
% to be implemented)
transNums      = NaN(3,3);
transNums(1,1) = length(strfind(stateStr, "CC"));   % Occurrences of Contraction to Contraction
transNums(1,2) = length(strfind(stateStr, "CS"));   % Occurrences of Contraction to Steady
transNums(1,3) = length(strfind(stateStr, "CE"));   % Occurrences of Contraction to Expansion
transNums(2,1) = length(strfind(stateStr, "SC"));   % Occurrences of Steady to Contraction
transNums(2,2) = length(strfind(stateStr, "SS"));   % Occurrences of Steady to Steady
transNums(2,3) = length(strfind(stateStr, "SE"));   % Occurrences of Steady to Expansion
transNums(3,1) = length(strfind(stateStr, "EC"));   % Occurrences of Expansion to Contraction
transNums(3,2) = length(strfind(stateStr, "ES"));   % Occurrences of Expansion to Steady
transNums(3,3) = length(strfind(stateStr, "EE"));   % Occurrences of Expansion to Expansion

% Find the number of times transitions occur from each state to any other
% state: the pairs in that state (leaving out the last pair), less those
% followed by a run start
numCSE2x    = zeros(3,1);
numCSE2x(1) = length(strfind(stateStr(1:end-1), "C")) - length(strfind(stateStr, "Cb"));
numCSE2x(2) = length(strfind(stateStr(1:end-1), "S")) - length(strfind(stateStr, "Sb"));
numCSE2x(3) = length(strfind(stateStr(1:end-1), "E")) - length(strfind(stateStr, "Eb"));

% Calculate the transition matrix (rows 2-4; row 1 stays NaN)
TM      = NaN(4,3);
TM(2,:) = transNums(1,:)./numCSE2x(1);
TM(3,:) = transNums(2,:)./numCSE2x(2);
TM(4,:) = transNums(3,:)./numCSE2x(3);

% TODO: Need to set up the Initial State row (TM row 1, left NaN above).
% The initial state row we do want to be influenced by a choice of weighting

%% Calculate the weighted TM (TMw)
% Set up weights choice: dx in cm for "depth", dt in kyr for "age", and 1
% per pair otherwise
if fitS.weighting == "depth"
    weights = nSRcountsArray(2,:);
elseif fitS.weighting == "age"
    weights = nSRcountsArray(3,:)./1000;   % Dividing by 1000 to convert units from yr to kyr
else
    weights = ones(size(nSRcountsArray(1,:)));
end

% Set up empty character and weight vectors, with one entry per stretch of
% consecutive pairs in the same state
blockStateStr = '';
blockWeights  = ones(0,0);

% Decide which category each nSR count belongs to (the same states as the
% TM loop above)
for i = 1:length(nSRs)
    if nSRs(i) >= 0.922 && nSRs(i) < 1.085
        pairState = 'S';
    elseif nSRs(i) >= 1.085 && nSRs(i) < inf
        pairState = 'C';
    elseif nSRs(i) >= 0 && nSRs(i) < 0.922
        pairState = 'E';
    else
        pairState = 'b';
    end

    % If the nSR is in the same category as the previous one (i.e., the
    % last one added to blockStateStr), just add the weights. If not, add a
    % new entry to blockStateStr and blockWeights
    if i ~= 1 && pairState == blockStateStr(end)
        blockWeights(end) = blockWeights(end) + weights(i);
    else
        blockStateStr = append(blockStateStr, pairState);
        blockWeights  = [blockWeights; weights(i)]; %#ok<AGROW>
    end
end

% Find the number of times for each transition: the self-transitions of a
% state are the summed weights of its stretches less one per stretch, and
% each change of state counts once
transNumsW      = NaN(3,3);
transNumsW(1,1) = sum(blockWeights(blockStateStr == 'C')) - sum(blockStateStr == 'C');
transNumsW(1,2) = length(strfind(blockStateStr, "CS"));   % Occurrences of Contraction to Steady
transNumsW(1,3) = length(strfind(blockStateStr, "CE"));   % Occurrences of Contraction to Expansion
transNumsW(2,1) = length(strfind(blockStateStr, "SC"));   % Occurrences of Steady to Contraction
transNumsW(2,2) = sum(blockWeights(blockStateStr == 'S')) - sum(blockStateStr == 'S');
transNumsW(2,3) = length(strfind(blockStateStr, "SE"));   % Occurrences of Steady to Expansion
transNumsW(3,1) = length(strfind(blockStateStr, "EC"));   % Occurrences of Expansion to Contraction
transNumsW(3,2) = length(strfind(blockStateStr, "ES"));   % Occurrences of Expansion to Steady
transNumsW(3,3) = sum(blockWeights(blockStateStr == 'E')) - sum(blockStateStr == 'E');

% Get the total (weighted) number of transitions from any given state
numCSE2xW = sum(transNumsW, 2);

% Calculate the transition matrix (weighted), with row 1 the fraction of
% the weighted transitions that start in each state
TMw      = NaN(4,3);
TMw(1,:) = numCSE2xW./sum(numCSE2xW);
TMw(2,:) = transNumsW(1,:)./numCSE2xW(1);
TMw(3,:) = transNumsW(2,:)./numCSE2xW(2);
TMw(4,:) = transNumsW(3,:)./numCSE2xW(3);

end
