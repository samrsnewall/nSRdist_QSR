function[binCounts] = makeWeightedBinCounts(data, weights, binEdges)
% This function makes a weighted histogram: for each bin set by binEdges it
% sums the weights of the data values that fall in that bin, so each value
% counts by its weight rather than by 1. Each bin includes its lower edge
% but not its upper edge, so a value equal to the last edge is not counted.
% binEdges are in the units of data (yr for the dt_yr bins of ARfitdists
% and IRfitdists); binCounts are in the units of weights (cm for dx_cm, yr
% for dt_yr, or a plain count for weights of one).
%
% S fields used: none

%% Sum the weights in each bin
% Initialise the bin index and the bin counts
binIdx    = 0;
binCounts = nan(1,length(binEdges)-1);

% For each bin...
for i = 1:length(binEdges)-1
    % Find the data values that fit in the bin
    inBinLog = (data >= binEdges(i) & data < binEdges(i+1));
    binIdx   = binIdx+1;   % Always equal to i

    % Find the count for that bin by summing the weights of those values
    % (dx_cm, dt_yr or ones, depending on the weighting)
    binCounts(binIdx) = sum(weights(inBinLog));
    inBinLog = [];   % Reset the logical vector
end
end
