function[sampData, numReplicates] = sampleWeightedReplicates(data, weights, weightRepDP, ...
    weightRepInflator, numSamples)
% This function draws numSamples values from the weighted-replicate dataset
% that makeWeightedReplicates would build, without building it. It is the
% fast alternative to makeWeightedReplicates followed by
% randsample(repData, numSamples, false), and draws from exactly the same
% distribution: the replicate counts are calculated in the same way, and the
% draw is uniform without replacement over the replicate positions. It is
% just done by drawing positions and looking up which unique value each
% position belongs to, rather than writing out every replicate first.
% numReplicates is the length the replicated dataset would have had.
% Units as in makeWeightedReplicates: in ARfitdists data are nSRs
% (unit-free, or SR in cm/kyr if fitS.non_normalized_SR) and weights are
% dx_cm or dt_yr.
%
% S fields used: none (ARfitdists passes its weightRepDataDP and
% weightRepInflator, and fitS.maxFitPoints as numSamples)

%% Check the inputs
% Warn if the data and weights are not the same length
if length(data) ~= length(weights)
    warning("Weights array not same size as data array - weights may not be properly aligned")
end

%% Find the replicate count of each unique value
% Round the data (to weightRepDP decimal places) and remove any values that
% rounded to zero, with their weights
data       = round(data, weightRepDP);
notZeroLog = data ~= 0;
data       = data(notZeroLog);
weights    = weights(notZeroLog);

% Accumulate the weights of each unique value, exactly as
% makeWeightedReplicates does, giving an integer replicate count per value
[uniqueData, ~, uniqueDataIdx] = unique(data);
uniqueWeights   = accumarray(uniqueDataIdx, weights);
uniqueRepCounts = round(uniqueWeights.*weightRepInflator);
numReplicates   = sum(uniqueRepCounts);

%% Draw the sample
% If the replicated dataset would be no bigger than the sample asked for,
% there is nothing to subsample: return the full replicated dataset
if numReplicates <= numSamples
    sampData = repelem(uniqueData, uniqueRepCounts);
    return
end

% Draw numSamples replicate positions (out of numReplicates) uniformly
% without replacement, as randsample would over the replicated vector
samplePositions = randperm(numReplicates, numSamples)';

% Find which unique value each position falls in: the replicated vector
% would hold value uniqueData(1) in positions 1 to uniqueRepCounts(1), then
% uniqueData(2) in the next uniqueRepCounts(2) positions, and so on
replicateEdges = [0; cumsum(uniqueRepCounts)];
valueIdx       = discretize(samplePositions, replicateEdges);
sampData       = uniqueData(valueIdx);
end
