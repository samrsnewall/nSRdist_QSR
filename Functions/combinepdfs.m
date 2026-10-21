function [master_invSRvals, master_invSRprobs] = combinepdfs(cells_invSRvals, cells_invSRprobs, lengths)
% combinepdfs  Combine per-scenario (or per-core) inverse-SR PDFs into a
%              single weighted PDF.
%
% Each scenario or core contributes a PDF defined on its own invSR grid.
% This function:
%   1. Discards scenarios/cores with empty PDFs (e.g. those that failed
%      reversal resolution).
%   2. Builds a master invSR grid spanning the union of all input grids,
%      using the grid spacing of the first non-empty PDF.
%   3. Places each input PDF at the correct position on the master grid,
%      scaled by the sediment length of that scenario/core.
%   4. Sums across all scaled PDFs and divides by the number of
%      scenarios/cores to produce the combined PDF.
%
% INPUTS
%   cells_invSRvals  - (cell array, s×1) Each cell holds the nSR value
%                      vector for one scenario or core. Cells may be
%                      empty for scenarios/cores that were rejected.
%   cells_invSRprobs - (cell array, s×1) Each cell holds the probability
%                      vector corresponding to cells_invSRvals.
%   lengths          - (numeric vector, s×1) Sediment length (cm) for each
%                      scenario or core, used as a weighting factor. NaN
%                      entries are dropped before use.
%
% OUTPUTS
%   master_invSRvals  - (numeric vector) Shared nSR axis for the combined PDF
%   master_invSRprobs - (numeric vector) Combined PDF probabilities on the
%                       master nSR axis
%
% See also: scenariopdfNorm, oneCoreScenarios

%% Remove scenarios/cores with empty PDFs
ind              = ~cellfun('isempty', cells_invSRvals);
cells_invSRvalsNE = cells_invSRvals(ind);
cells_invSRprobs  = cells_invSRprobs(ind);
numscenarios      = length(cells_invSRvalsNE);

if numscenarios == 0
    error("None of the cores used passed all the criteria, hence there are no pdfs of invSR to combine")
end

%% Remove NaN lengths
lengths = lengths(~isnan(lengths));

%% Build a master nSR grid spanning all input PDFs
% Grid spacing is taken from the first non-empty PDF; all inputs are
% assumed to share the same spacing (set in scenariopdfNorm).
min_invSRval = min(cellfun(@min, cells_invSRvalsNE));
max_invSRval = max(cellfun(@max, cells_invSRvalsNE));
interval     = cells_invSRvalsNE{1}(2) - cells_invSRvalsNE{1}(1);
master_invSRvals = (min_invSRval:interval:(max_invSRval + interval/2))';

%% Place each PDF onto the master grid and scale by sediment length
m_invSRprobs = zeros(length(master_invSRvals), numscenarios);
for i = 1:numscenarios
    % Find the index on the master grid where this PDF begins.
    idx = knnsearch(master_invSRvals, cells_invSRvalsNE{i}(1));
    m_invSRprobs(idx:idx + length(cells_invSRvalsNE{i}(:)) - 1, i) = ...
        cells_invSRprobs{i}(:) .* lengths(i);
end

%% Average across all scenarios/cores
master_invSRprobs = sum(m_invSRprobs, 2) ./ numscenarios;

end
