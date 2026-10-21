function[invGamPDF, GamPDF, invGamFitInfo] = fitInvGamma(weightedC, x, numCpairs)
% This function fits an inverse gamma distribution (invGam) to the data
% weightedC by fitting a gamma distribution to 1./weightedC by maximum
% likelihood (gamfit). It gives the inverse gamma pdf on the grid x
% (invGamPDF, transformed from the gamma pdf on 1./x with px_to_pfx) and the
% gamma pdf on 1./x (GamPDF, for checking). weightedC can be weighted
% replicates (see makeWeightedReplicates), and numCpairs is the number of
% observations before replication.
%
% invGamFitInfo holds the shape (alpha) and scale (beta) of the gamma fitted
% to 1./weightedC, the number of parameters, the NLL of weightedC (with the
% Jacobian correction 2*sum(log(weightedC)) added so it is the NLL in x
% space), the AIC, the BIC, and BICtaeheefix, the BIC with the NLL divided
% by the replication factor length(weightedC)./numCpairs and numCpairs as
% the number of observations. weightedC and x are nSR (unit-free; SR in
% cm/kyr if the callers' fitS.non_normalized_SR), so the gamma is fitted to
% 1/nSR.
%
% fitS fields used: none

%% Fit a gamma distribution to the inverse of the data
% Set the number of iterations (default=100 often doesn't converge)
fitOptions = statset('Display', 'off', 'MaxIter', 200, 'TolFun', 1e-6);

% Ensure the data is in column vector form
if ~iscolumn(weightedC)
    weightedC = weightedC';
end

% Take the inverse of the data, fit the gamma distribution to it, and get
% its shape and scale
invData   = 1./weightedC;
gamParams = gamfit(invData, "Options", fitOptions);
alpha     = gamParams(1);
beta      = gamParams(2);

%% Create the gamma pdf on the inverse of the grid x
invx   = 1./x;
GamPDF = gampdf(invx, alpha, beta)';

%% Find the likelihood of the fitted distribution
% Find the NLL of the inverse data under the fitted gamma, and add the
% Jacobian correction 2*sum(log(weightedC)) to give the NLL of weightedC
NLLinvData   = gamlike([alpha, beta], invData);
jacobianTerm = 2*sum(log(weightedC));
NLL          = NLLinvData + jacobianTerm;

%% Find the inverse gamma pdf on the grid x
% Transform the gamma pdf on invx to x space (px_to_pfx gives it on the grid
% sort(1./invx))
invFH               = @(x) 1./x;
[invGamX, invGamPx] = px_to_pfx(invx, GamPDF, invFH);

% Interpolate back to the desired x vector spacing (both grids rounded to 10
% decimal places, so that the grid that comes back from 1./(1./x) matches x)
invGamXRounded = round(invGamX, 10);
xRounded       = round(x, 10);
invGamPDF      = interp1(invGamXRounded, invGamPx, xRounded)';

%% Correct the BIC for the replication
% Calculate NLL and correct for number of repetitions, as suggested by
% Taehee (divide NLL by number of repetitions), and find the BIC with
% numCpairs observations (2 parameters)
replicationFactor = length(weightedC)./numCpairs;
NLLtaeheeFix      = NLL./replicationFactor;
BICtaeheeFix      = 2*NLLtaeheeFix + 2*log(numCpairs);

%% Set up output structure with important values
% AIC = 2k - 2ln(L), therefore = 2k + 2*NLL
invGamFitInfo.alpha                 = alpha;
invGamFitInfo.beta                  = beta;
invGamFitInfo.NumParams             = 2;
invGamFitInfo.NegativeLogLikelihood = NLL;
invGamFitInfo.AIC                   = 2*invGamFitInfo.NumParams + 2*NLL;
invGamFitInfo.BIC                   = log(length(weightedC))*invGamFitInfo.NumParams + 2*NLL;
invGamFitInfo.BICtaeheefix          = BICtaeheeFix;

%% Check on likelihood
% DEBUG: work out the NLL again straight from the inverse gamma pdf, as a
% check on NLL above (NLLcheck should equal NLL; none of this is output)
GamPDFcheckFH    = @(y) gampdf(y, alpha, beta);
invGamPDFcheckFH = @(x) gampdf(1./x, alpha, beta) .* (1./x.^2);
x                = weightedC;
NLLcheck         = -sum(log(invGamPDFcheckFH(x)));
end
