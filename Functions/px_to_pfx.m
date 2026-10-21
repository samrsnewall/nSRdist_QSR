function[fx, pfx] = px_to_pfx(x, px, f)
% This function transforms a pdf p(x), given on the grid x, into the pdf of
% f(x) for a monotonic function f (e.g. @log, or @(x) 1./x), giving the
% densities pfx on the grid fx = sort(f(x)). Each density is divided by the
% width its bin takes up after applying f, the result is interpolated onto
% fx, and the pdf is normalised so that trapz(fx, pfx) = 1. Grid points
% below 1e-5 are left out of the transformation (the fx points outside the
% transformed range get the density of the nearest point inside it), and if
% x is not equally spaced, p(x) is first interpolated onto an equally spaced
% grid. fx is in the units of f(x): with x in nSR (unit-free) and f = @log,
% fx is log(nSR).
%
% S fields used: none

%% Keep the grid points above zero
% Keep x and px where x is at least 1e-5
xPositive  = x(x>=1e-5);
pxPositive = px(x>=1e-5);

%% Test whether the x is on equal spacing
% Find the unique values of dx, how different each is from the first, and
% how many of them are actually significant (more than 1e-5 apart,
% counting the first)
uniqueDx         = unique(diff(x));
uniqueDxResidual = abs(uniqueDx - uniqueDx(1));
distinctDxLog    = uniqueDxResidual(1:end)>1e-5;
distinctDxLog(1) = true;
numDistinctDx    = sum(distinctDxLog);

% If x is not equally spaced, convert it to equal spacing (bin width 5e-4)
% by interpolating the pdf; otherwise the bin width is the spacing of x
if numDistinctDx ~= 1
    xBinWidth     = 5e-4;
    xEqualSpaced  = min(xPositive):xBinWidth:max(xPositive);
    pxEqualSpaced = interp1(xPositive, pxPositive, xEqualSpaced);

    % Assign these to the variables used in the next part of the code
    xPositive  = xEqualSpaced;
    pxPositive = pxEqualSpaced;
else
    % FIXME: in this branch xPositive and pxPositive keep the orientations
    % of x and px, so an x and px of different orientations (one a row, one
    % a column) give a square matrix for pfxBins below instead of a vector
    % (the branch above always gives rows)
    xBinWidth = uniqueDx(1);
end

%% Apply correction for differing bin sizes after applying function
% Find the bin edges and the bin width in fx space
% FIXME: grid points between 1e-5 and xBinWidth./2 make xPositive -
% xBinWidth./2 negative, so with f = @log fxHalfStepDown is complex and
% interp1 below errors
fxHalfStepDown = f(xPositive - xBinWidth./2);
fxHalfStepUp   = f(xPositive + xBinWidth./2);
fxBinSizes     = abs(fxHalfStepDown - fxHalfStepUp);

% Divide the prob by the bin size in fx space
pfxBins = pxPositive./fxBinSizes;

% Find the fx values that apply to the bins (the middle of each bin in fx
% space)
fxBins = (fxHalfStepDown + fxHalfStepUp)./2;

%% Interpolate onto the original grid after applying the function
% Create the new spacing from the original x, and find the points outside
% the range of fxBins (interp1 gives NaN there, so they are filled in below)
fxGrid      = sort(f(x));
tooBigLog   = fxGrid > max(fxBins);
tooSmallLog = fxGrid < min(fxBins);

% Interpolate to the new spacing, and give the points outside the range of
% fxBins the density of the nearest point inside it
pfxGrid              = interp1(fxBins, pfxBins, fxGrid);
pfxGrid(tooBigLog)   = pfxGrid(find(tooBigLog, 1, 'first') -1);
pfxGrid(tooSmallLog) = pfxGrid(find(tooSmallLog, 1, 'last') +1);

%% Normalise the prob vector
% Divide by the area under the curve so it integrates to 1 over fxGrid
AUC     = trapz(fxGrid, pfxGrid);
pfxNorm = pfxGrid./AUC;

% Assign to values to be output
fx  = fxGrid;
pfx = pfxNorm;
end
