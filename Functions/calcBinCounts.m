function[obsCD, expCountsArr] = calcBinCounts(data, binEdges, desiredSum, pdfVs)
% calcBinCounts  Compute observed and expected bin counts for a chi-squared test.
%
% Bins the observed data using histcounts, then downscales the counts by
% the ratio numel(data)/desiredSum so that the total observed count equals
% desiredSum. This rescaling is used when data is a set of weighted
% replicates (whose total may be much larger than the true sample size) and
% desiredSum is the true number of unique observations; it prevents the
% chi-squared test from being over-powered by replicated data.
%
% For each candidate PDF in pdfVs, the expected bin count is computed from
% the PDF's CDF: the probability mass within each bin is obtained from
% differences of CDF values at the bin edges, then multiplied by
% desiredSum to give an expected count.
%
% INPUTS
%   data        - (numeric vector) Observed data values (may be weighted
%                 replicates from makeWeightedReplicates)
%   binEdges    - (numeric vector, length B+1) Edges of the B histogram bins
%   desiredSum  - (scalar) Target total count; observed counts are
%                 downscaled to this value, and expected counts are computed
%                 to sum to this value
%   pdfVs       - (cell array, P×1) Candidate PDFs; each cell is a struct
%                 with fields:
%                   .x      (numeric vector) PDF evaluation points
%                   .px     (numeric vector) Probability density values
%                   .cdf_x  (numeric vector) Normalised CDF at each .x point
%                           (must be pre-computed by the caller)
%
% OUTPUTS
%   obsCD        - (1 × B numeric vector) Observed bin counts, downscaled
%                  and rounded to integers
%   expCountsArr - (P × B numeric matrix) Expected bin counts for each PDF
%
% See also: chi2_dataVStwopdfVECs, makeWeightedReplicates

%Find how many observations are in each bin
obsC  = histcounts(data, binEdges);
divisor = numel(data)/desiredSum;
obsCD = round(obsC./divisor);

%Get expected counts for each pdf
expCountsArr = NaN(size(pdfVs, 1), length(binEdges)-1);
for i = 1:size(pdfVs,1)
    pdfi = pdfVs{i}; %Choose pdf
    cdfAtEdges = interp1(pdfi.x, pdfi.cdf_x, binEdges(2:end-1)); %Get cdf values at each internal binEdge
    expProbsi = diff([0, cdfAtEdges, 1]); %Find difference of cdf values at each binEdge (bottom must be 0, top must be 1), this provides what fraction of probability is within this bin
    expCountsi = expProbsi.*desiredSum; % Multiply this by desired sum
    expCountsArr(i,:) = expCountsi; %Save results to array
end
end