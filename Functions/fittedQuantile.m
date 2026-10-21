function q = fittedQuantile(fitInfo, family, p)
% fittedQuantile  Quantile function (inverse CDF) of a fitted nSR distribution.
%
% Closed form for LN, Gam and invGam. The two-component MLN has no closed-form
% inverse, so its CDF is evaluated on a dense log(nSR) grid and inverted by
% interpolation (grid spacing 0.001 in log nSR, i.e. ~0.1 % in nSR, which is
% far below plotting resolution).
%
% INPUTS
%   fitInfo - fit-information struct, e.g. d.S1.RSR500AR.MLN.fitInfo
%   family  - "MLN", "LN", "Gam" or "invGam"
%   p       - probabilities in (0, 1)
%
% OUTPUT
%   q       - nSR quantiles (linear space), same size as p
%
% See also: fittedCDF, qqPairs

family = string(family);
p = double(p);
if any(p(:) <= 0 | p(:) >= 1)
    error("fittedQuantile: probabilities must lie strictly inside (0, 1)")
end

switch family
    case "LN"
        % Single log-normal: log X ~ N(mu, sigma^2)
        mu  = fitInfo.mu(1);
        sig = sqrt(squeeze(fitInfo.Sigma(1)));
        q   = exp(mu + sig .* norminv(p));

    case "Gam"
        q = gaminv(p, fitInfo.alpha, fitInfo.beta);

    case "invGam"
        % X = 1/Y, Y ~ Gamma(alpha, beta):  P(X <= q) = p  <=>  P(Y >= 1/q) = p
        %                                              <=>  gamcdf(1/q) = 1 - p
        q = 1 ./ gaminv(1 - p, fitInfo.alpha, fitInfo.beta);

    case "MLN"
        % Grid inversion. Extend the grid until it covers the requested tail
        % probabilities, so extreme p values are never extrapolated.
        lgrid = -12:0.001:12;
        Fgrid = fittedCDF(fitInfo, "MLN", exp(lgrid));
        % Make strictly increasing for interp1 (drop flat tails at 0 and 1)
        keep  = [true, diff(Fgrid) > 0];
        Fgrid = Fgrid(keep);
        lgrid = lgrid(keep);
        if min(p(:)) < Fgrid(1) || max(p(:)) > Fgrid(end)
            error("fittedQuantile: MLN grid does not cover requested p range (%.3g to %.3g)", min(p(:)), max(p(:)))
        end
        q = exp(interp1(Fgrid, lgrid, p, "linear"));

    otherwise
        error("fittedQuantile: unknown family '%s' (expected MLN, LN, Gam or invGam)", family)
end

end
