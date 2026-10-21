function[mu, Var] = muVarPDFVec(pdfVEC)
% This function calculates the mean and variance of a pdf given on a grid
% (pdfVEC.x, with densities pdfVEC.px) by numerical integration by summing
% areas (trapz) of x.*px and x.^2.*px. The pdf is not normalised here, so
% the results only mean something if it integrates to about 1. mu and Var
% are in the units of pdfVEC.x: nSR is unit-free, and for a pdf on log(nSR)
% (the .lnSR structs) they are the mean and variance of log(nSR).
%
% S fields used: none

%% Calculate the mean and variance
% Get the grid and the densities on it
x  = pdfVEC.x;
px = pdfVEC.px;

% Integrate x*p(x) for the mean, and x^2*p(x) minus the squared mean for
% the variance
mu  = trapz(x, x.*px);
Var = trapz(x, x.*x.*px) - mu^2;
end
