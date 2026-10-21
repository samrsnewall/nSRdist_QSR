function[repData] = makeWeightedReplicates(data, weights, weightRepDP, weightRepInflator)
% This function turns weighted data into an unweighted vector of replicates
% that approximates the same weighting, because the fitters (fitgmdist,
% gamfit etc) need individual observations and cannot take weights or bin
% counts. The data are rounded to weightRepDP decimal places, the weights of
% each unique rounded value are summed, multiplied by weightRepInflator and
% rounded to an integer, and each unique value is repeated that many times.
% In ARfitdists and IRfitdists data are nSRs (unit-free, or SR in cm/kyr if
% fitS.non_normalized_SR) and weights are dx_cm or dt_yr, so
% weightRepInflator is in replicates per cm or per yr.
%
% S fields used: none (ARfitdists and IRfitdists pass
% fitS.PooledRuns.weightRepDP and weightRepInflator, or the fitS.OneRun ones)

%% Check the inputs
% Warn if the data and weights are not the same length
if length(data) ~= length(weights)
    warning("Weights array not same size as data array - weights may not be properly aligned")
end

%% Round the data and remove zeros
% Round the data (to weightRepDP decimal places) so that there are fewer
% unique values, to help count the weightings
data = round(data, weightRepDP);

% Remove any values that rounded down to zero (data should be strictly
% positive), and their weights
notZeroLog = data ~= 0;
data       = data(notZeroLog);
weights    = weights(notZeroLog);

%% Replicate each unique value by its weight
% Find the unique values of the data, with the index of the unique value
% each data value belongs to
[uniqueData, ~, uniqueDataIdx] = unique(data);

% Combine the weights of each unique value
uniqueWeights = accumarray(uniqueDataIdx, weights);

% Multiply the combined weight of each unique value by weightRepInflator
% (useful to avoid rounding to 0) and then round it to an integer. A value
% whose weight still rounds to 0 is dropped; increase weightRepInflator to
% reduce this
uniqueRepCounts = round(uniqueWeights.*weightRepInflator);

% Replicate the unique values according to their integer weighting
repData = repelem(uniqueData, uniqueRepCounts);

% DEBUG: Test to see if the approach above makes much of a difference to
% the very simple approach of
% weight_R = round(weights.*weightRepInflator);
% repData2 = repelem(data, weight_R);
%
% binEdges = 0:0.05:10;
%
% figure;
% subplot(2,1,1)
% histogram(repData, 'BinEdges', binEdges)
% title("repData")
% subplot(2,1,2)
% histogram(repData2, 'BinEdges', binEdges)
% title("repData2")
%
% if any(uniqueRepCounts == 0) || any(weight_R == 0)
%     lostWeight1 = sum(uniqueWeights(uniqueRepCounts == 0));
%     lostWeight2 = sum(weights(weight_R == 0));
%     warning("repData lost " + num2str(lostWeight1) + " weight ; repData2 lost " + num2str(lostWeight2) + " weight")
% end

end
