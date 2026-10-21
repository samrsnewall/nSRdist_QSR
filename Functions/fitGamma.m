function[GamPDF, GamFitInfo] = fitGamma(weightedC, x, numCpairs)
% This function fits a gamma distribution (Gam) to the data weightedC by
% maximum likelihood (gamfit) and gives its pdf on the grid x. weightedC can
% be weighted replicates (see makeWeightedReplicates), and numCpairs is the
% number of observations before replication.
%
% GamFitInfo holds the shape (alpha) and scale (beta), the number of
% parameters, the NLL, the AIC, the BIC, and BICtaeheefix, the BIC with the
% NLL divided by the replication factor length(weightedC)./numCpairs and
% numCpairs as the number of observations. weightedC and x are nSR
% (unit-free; SR in cm/kyr if the callers' fitS.non_normalized_SR).
%
% fitS fields used: none

%% Fit a gamma distribution to the data
% Set the number of iterations (default=100 often doesn't converge)
fitOptions = statset('Display', 'off', 'MaxIter', 200, 'TolFun', 1e-6);

% Ensure the data is in column vector form
if ~iscolumn(weightedC)
    weightedC = weightedC';
end

% Fit the gamma distribution, and get its shape and scale
gamParams = gamfit(weightedC, "Options", fitOptions);
alpha     = gamParams(1);
beta      = gamParams(2);

%% Create the gamma pdf on the grid x
GamPDF = gampdf(x, alpha, beta)';

%% Find the likelihood of the fitted distribution
% Find the NLL of the data under the fitted gamma
NLL = gamlike([alpha, beta], weightedC);

% Calculate NLL and correct for number of repetitions, as suggested by
% Taehee (divide NLL by number of repetitions), and find the BIC with
% numCpairs observations (2 parameters)
replicationFactor = length(weightedC)./numCpairs;
NLLtaeheeFix      = NLL./replicationFactor;
BICtaeheeFix      = 2*NLLtaeheeFix + 2*log(numCpairs);

%% Set up output structure with important values
% AIC = 2k - 2ln(L), therefore = 2k + 2*NLL
GamFitInfo.alpha                 = alpha;
GamFitInfo.beta                  = beta;
GamFitInfo.NumParams             = 2;
GamFitInfo.NegativeLogLikelihood = NLL;
GamFitInfo.AIC                   = 2*GamFitInfo.NumParams + 2*NLL;
GamFitInfo.BIC                   = log(length(weightedC))*GamFitInfo.NumParams + 2*NLL;
GamFitInfo.BICtaeheefix          = BICtaeheeFix;
end
