function qq = qqPairs(methodS, family, opts)
% qqPairs  Empirical-vs-fitted quantile pairs for a Q-Q plot of one nSR
%          estimation method against one fitted distribution family.
%
% The empirical quantiles come from methodS.weightedC, the replicated
% (depth- or age-weighted) data vector that was actually passed to the
% fitting routines, so they are the weighted empirical quantiles.
%
% Plotting positions are set by the EFFECTIVE sample size, not by the length
% of weightedC: p_i = (i - 0.5) / n_eff with n_eff = round(methodS.numCpairs),
% the same n that BICeff uses. The most extreme quantile examined is therefore
% p = 0.5/n_eff, i.e. as far into the tails as ~n_eff independent observations
% can support. Using length(weightedC) instead would push the plot into
% percentiles that the replication has manufactured rather than observed.
%
% Optionally returns a pointwise 95 % envelope for the empirical quantile
% under the fitted model, from the Beta(i, n_eff - i + 1) distribution of the
% i-th order statistic's probability level. This assumes n_eff independent
% draws - the assumption Reviewer 1 questioned - so label it as approximate.
%
% INPUTS
%   methodS - method result struct, e.g. d.S1.RSR500AR (needs .weightedC,
%             .numCpairs and .<family>.fitInfo)
%   family  - "MLN", "LN", "Gam" or "invGam"
%   opts    - (optional struct)
%             .nEff       override effective sample size (default round(numCpairs))
%             .envelope   logical, compute 95 % envelope (default true)
%             .envLevel   two-sided level for the envelope (default 0.95)
%
% OUTPUT (struct)
%   qq.p      plotting positions (n_eff x 1)
%   qq.qEmp   empirical nSR quantiles at p
%   qq.qFit   fitted nSR quantiles at p
%   qq.lo,.hi envelope bounds on the empirical quantile (if requested)
%   qq.nEff   effective sample size used
%   qq.nRep   length of the replicated data vector (for reference)
%   qq.family, qq.BICeff
%
% WHY NOT MATLAB'S qqplot / probplot
% MATLAB's native Q-Q plot functions don't easily support some of our needs,
% so we write a bespoke function here. In particular, they don't support:
%   * drawing the plots against ln(nSR) instead of nSR (there is no built-in
%     "log gamma" or "log inverse gamma" distribution, and the MLN
%     gmdistribution has no icdf)
%   * setting an effective sample size different from the vector length
%     (qqplot takes n = length(x), i.e. the replicated weightedC)
%
% See also: fittedCDF, fittedQuantile, qq_tailGOF

if nargin < 3, opts = struct(); end
if ~isfield(opts, "envelope"), opts.envelope = true;  end
if ~isfield(opts, "envLevel"), opts.envLevel = 0.95;  end

family = string(family);

% Replicated data vector (older fit files stored it in a cell)
wc = methodS.weightedC;
if iscell(wc), wc = cell2mat(wc(:)); end
wc = wc(:);
wc = wc(isfinite(wc) & wc > 0);

if isfield(opts, "nEff") && ~isempty(opts.nEff)
    nEff = round(opts.nEff);
else
    nEff = round(mean(methodS.numCpairs));
end

% Plotting positions from the effective sample size
i = (1:nEff)';
p = (i - 0.5) ./ nEff;

qq.p      = p;
qq.qEmp   = quantile(wc, p);
qq.qEmp   = qq.qEmp(:);
qq.qFit   = fittedQuantile(methodS.(family).fitInfo, family, p);
qq.qFit   = qq.qFit(:);
qq.nEff   = nEff;
qq.nRep   = numel(wc);
qq.family = family;
qq.BICeff = methodS.(family).fitInfo.BICtaeheefix;

if opts.envelope
    a  = (1 - opts.envLevel) / 2;
    pl = betainv(a,     i, nEff - i + 1);
    ph = betainv(1 - a, i, nEff - i + 1);
    qq.lo = fittedQuantile(methodS.(family).fitInfo, family, pl);
    qq.hi = fittedQuantile(methodS.(family).fitInfo, family, ph);
    qq.lo = qq.lo(:);
    qq.hi = qq.hi(:);
    qq.envLevel = opts.envLevel;
end

end
