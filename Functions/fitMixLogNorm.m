function[mixLogNormPDF, mixNormPDF, MLNfitInfo] = fitMixLogNorm(weightedC, x, numComponents, ...
    numRandomStarts, numCpairs)
% This function fits a mixed log-normal (MLN) of numComponents components to
% the data weightedC, by fitting a Gaussian mixture to log(weightedC) with
% fitgmdist (numRandomStarts random starts, fitgmdist's 'Replicates') and
% turning each Gaussian component into the log-normal with the same mu and
% sigma. numComponents = 1 gives a single log-normal (LN). weightedC can be
% weighted replicates (see makeWeightedReplicates), and numCpairs is the
% number of observations before replication.
%
% mixLogNormPDF is the MLN pdf on the grid x (columns x and density), and
% mixNormPDF is the Gaussian mixture pdf on log(x) (columns log(x) and
% density), for checking the fit in log space. MLNfitInfo holds the number
% of parameters, the fitted mu, Sigma (variances) and ComponentProportion,
% the NLL of weightedC (NegativeLogLikelihood, with the Jacobian correction
% sum(log(weightedC)) added so it is the NLL in linear space), the BIC, and
% BICtaeheefix, the BIC with the NLL divided by the replication factor
% length(weightedC)./numCpairs and numCpairs as the number of observations.
% weightedC and x are nSR (unit-free; SR in cm/kyr if the callers'
% fitS.non_normalized_SR); mu, Sigma and the grid of mixNormPDF are in
% log(nSR).
%
% fitS fields used: none (the callers pass fitS.mlnReps or
% fitS.OneRun.MLNReps as numRandomStarts)

%% Fit a mixed normal to the log of the data
% Set the number of iterations (default=500 often doesn't converge)
fitOptions = statset('Display', 'off', 'MaxIter', 1000, 'TolFun', 1e-6);

% Ensure the data is in column vector form
if ~iscolumn(weightedC)
    weightedC = weightedC';
end

% Take the logarithm of the data and fit a mixed normal to it
logData = log(weightedC);
gmfit   = fitgmdist(logData, numComponents, "Options", fitOptions, ...
    'Replicates', numRandomStarts, 'Start','randSample');

%% Find the pdf of each component
% Initialise the std deviation vector, and the pdf of each component as a
% log-normal on x and as a normal on log(x)
sigma                = NaN(1,numComponents);
componentLogNormPDFs = zeros(length(x), numComponents);
logx                 = log(x);
componentNormPDFs    = zeros(length(logx), numComponents);

% Pull out the std deviation of each component Gaussian (the square root of
% its variance), and find the log-normal pdf on x and the normal pdf on
% log(x) using the mean and std deviation of the component Gaussian as the
% parameters
mu = gmfit.mu;
for i = 1:numComponents
    sigma(i)                  = sqrt(gmfit.Sigma(:,:,i));
    componentLogNormPDFs(:,i) = lognpdf(x, mu(i), sigma(i));
    componentNormPDFs(:,i)    = normpdf(logx, mu(i), sigma(i));
end

%% Find the mixed log-normal
% Sum the component log-normals after multiplying each by its mixing
% proportion
componentProportion = gmfit.ComponentProportion;
mixLogNormPx        = sum(componentLogNormPDFs.*componentProportion, 2);

% Flag an MLN that peaks at the first grid point
if find(mixLogNormPx == max(mixLogNormPx)) == 1
    disp('the max value is the first index... This is suspicious, double check if this is meant to be the case')
end

% Sum the component normals on log(x) in the same way, and put each pdf
% next to its grid
mixNormPx     = sum(componentNormPDFs.*componentProportion, 2);
mixLogNormPDF = [x', mixLogNormPx];
mixNormPDF    = [logx', mixNormPx];

%% Calculate the negative log likelihood in linear space
% Count the parameters (mu and sigma of the first component, plus mu,
% sigma and a mixing proportion for each further one)
numParams = 2+(numComponents-1)*3;

% Add sum(log(weightedC)) to the NLL of the fit to log(weightedC) (the
% Jacobian correction to convert back to linear space), and find the BIC
jacobianTerm = sum(log(weightedC));
NLL          = gmfit.NegativeLogLikelihood + jacobianTerm;
BIC          = 2*NLL + numParams*log(length(weightedC));

% Calculate NLL and correct for number of repetitions, as suggested by
% Taehee (divide NLL by number of repetitions), and find the BIC with
% numCpairs observations
replicationFactor = length(weightedC)./numCpairs;
NLLtaeheeFix      = NLL./replicationFactor;
BICtaeheeFix      = 2*NLLtaeheeFix + numParams*log(numCpairs);

%% Set up output structure
MLNfitInfo.NumParams             = numParams;
MLNfitInfo.mu                    = gmfit.mu;
MLNfitInfo.Sigma                 = gmfit.Sigma;
MLNfitInfo.ComponentProportion   = gmfit.ComponentProportion;
MLNfitInfo.NegativeLogLikelihood = NLL;
MLNfitInfo.BIC                   = BIC;
MLNfitInfo.BICtaeheefix          = BICtaeheeFix;
end
