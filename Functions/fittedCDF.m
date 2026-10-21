function F = fittedCDF(fitInfo, family, x)
% fittedCDF  Closed-form CDF of a fitted nSR distribution, evaluated at x.
%
% Uses the parameters stored in <method>.<family>.fitInfo by ARfitdists
% (via fitMixLogNorm, fitGamma, fitInvGamma), so the CDF is exactly that
% of the distribution whose PDF is stored on the d.x grid - no numerical
% integration of the PDF grid is involved.
%
% INPUTS
%   fitInfo - fit-information struct, e.g. d.S1.RSR500AR.LN.fitInfo
%   family  - "MLN", "LN", "Gam" or "invGam"
%   x       - nSR values (linear space, > 0) at which to evaluate the CDF
%
% OUTPUT
%   F       - CDF values, same size as x
%
% PARAMETERISATION (matches the fitting functions)
%   MLN / LN : fitgmdist on log(nSR); fitInfo.mu (K x 1), fitInfo.Sigma
%              (1 x 1 x K, variances), fitInfo.ComponentProportion (1 x K).
%              F(x) = sum_k pi_k * Phi( (log x - mu_k) / sigma_k )
%   Gam      : gamfit on nSR; fitInfo.alpha (shape), fitInfo.beta (scale),
%              so F(x) = gamcdf(x, alpha, beta)
%   invGam   : gamfit on 1/nSR; fitInfo.alpha, fitInfo.beta describe
%              Y = 1/X ~ Gamma(alpha, beta). Then
%              F(x) = P(X <= x) = P(Y >= 1/x) = 1 - gamcdf(1/x, alpha, beta)
%
% See also: fittedQuantile, qqPairs, ARfitdists

family = string(family);
x = double(x);
F = NaN(size(x));
pos = x > 0;

switch family
    case {"MLN", "LN"}
        mu    = fitInfo.mu(:);                                % K x 1
        sig   = sqrt(squeeze(fitInfo.Sigma));                 % K x 1 (std devs)
        sig   = sig(:);
        w     = fitInfo.ComponentProportion(:);               % K x 1
        lx    = log(x(pos));
        Fpos  = zeros(size(lx));
        for k = 1:numel(w)
            Fpos = Fpos + w(k) .* normcdf((lx - mu(k)) ./ sig(k));
        end
        F(pos)  = Fpos;
        F(~pos) = 0;

    case "Gam"
        F(pos)  = gamcdf(x(pos), fitInfo.alpha, fitInfo.beta);
        F(~pos) = 0;

    case "invGam"
        F(pos)  = 1 - gamcdf(1 ./ x(pos), fitInfo.alpha, fitInfo.beta);
        F(~pos) = 0;

    otherwise
        error("fittedCDF: unknown family '%s' (expected MLN, LN, Gam or invGam)", family)
end

end
