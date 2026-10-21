function[outS] = IRfitdists(nSRcounts, coreLog, numruns, x, applyLin2014AgeFilter, fitS)
% This function fits the same four distributions as ARfitdists (a
% 2-component mixed log-normal (MLN), a log-normal (LN), a gamma (Gam) and an
% inverse gamma (invGam)) to one randomly chosen run from each chosen core
% (coreLog), numruns times, giving a distribution of fit parameters that
% shows how much the fits vary from run to run, where ARfitdists gives one
% fit to all runs pooled. Each iteration's data are cleaned and weighted as
% in ARfitdists (fitS.weighting), and if applyLin2014AgeFilter is true only
% the pairs whose dt_yr is at least the lower value of fitS.Lin2014AgeFilter
% and below the upper value are kept, as in ARfitdists (fitData does this
% for BSampLinIR). nSR is unit-free (SR is in cm/kyr if
% fitS.non_normalized_SR); dx in cm; dt in yr.
%
% nSRcounts is one dataT column (one 3-row nSR matrix per core: nSR, dx_cm,
% dt_yr, each run starting with a NaN column); chosen cores with an empty
% cell are dropped, and every other chosen core needs at least numruns runs.
% x is the nSR grid the pdfs are given on. outS holds, for each iteration
% (one column or cell each), the runs drawn before any cleaning, each
% starting with its NaN run-start column
% (.OneRunDatas), the data the fitters saw (.weightedC), the number of nSR
% values (.numCpairs, which the fitters use as the effective number of
% observations), the weighted bin counts of dt_yr (.agediffsWbC), the summed
% dx (cm) and dt (yr) (.sedLength, .sedTimeSpan), and for each of .MLN,
% .invGam, .Gam and .LN the fit results (.fits) and the .nSR and .lnSR pdfs
% with their mean and variance. It also holds one chi-squared table per
% distribution (.chiStats, NaN unless fitS.run_chi2gof) and
% .Lin2014AgeFiltering.
%
% fitS fields used: useParallelComp, non_normalized_SR,
% Lin2014AgeFilter, OneRun.weightRepDP, OneRun.weightRepInflator, weighting,
% resampleData, fitDists, OneRun.MLNReps, run_chi2gof (and those
% chi2_dataVStwopdfVECs reads, with dispChi2 set here)

%% Set up
% Set up the log(nSR) grid, the number of grid points and of pdfs, and the
% dt_yr bin edges (in yr)
logx          = sort(log(x));
numX          = numel(x);
numPdfs       = 4;
dtBinEdges_yr = 0:100:10000;

% Remove the chosen cores that returned empty nSR matrices, and initialise
% which run is drawn from each core in each iteration
nSRcountsChosen = nSRcounts(coreLog);
emptyCoreLog    = cellfun(@isempty, nSRcountsChosen);
nSRcountsClean  = nSRcountsChosen(~emptyCoreLog);
chosenRunIdx    = NaN(numruns, length(nSRcountsClean));

%% Sample one run per core, for each of numruns iterations
% Put one run from each core side by side for each iteration, each run
% starting with its NaN column (how ARfitdists's pooled array looks)
OneRunDatas = cell(1, numruns);
for i = 1:numruns
    OneRunData = [];
    for j = 1:length(nSRcountsClean)
        core_nSRcounts = nSRcountsClean{j};

        % Find the NaN run-start columns, and the last column of each run
        runStartIdx = find(isnan(core_nSRcounts(1,:)));
        runLastIdx  = [runStartIdx(2:end) - 1, size(core_nSRcounts, 2)];

        % Randomise which run is drawn from this core in each iteration (done
        % once, in the first iteration, without replacement, so each
        % iteration draws a different run)
        if i == 1
            chosenRunIdx(:,j) = randperm(length(runStartIdx), numruns);
        end

        % Get the chosen run of this core, and flag it if it has zero or
        % negative nSRs (these are removed before fitting)
        runIdx        = chosenRunIdx(i,j);
        run_nSRcounts = core_nSRcounts(1:3, runStartIdx(runIdx):runLastIdx(runIdx));
        if sum(run_nSRcounts(1,:) <= 0) ~= 0
            disp("check this core")
        end

        % Add it to this iteration's data
        OneRunData = cat(2, OneRunData, run_nSRcounts);
    end
    OneRunDatas{i} = OneRunData;
end

%% Set up the parallel pool
% Open a pool of 4 workers if fitS.useParallelComp and none is open;
% otherwise stop MATLAB opening one automatically (restored at the end)
% FIXME: the fitting loop below is a plain for loop (no parfor) and the
% fitters do not ask for a pool, so fitS.useParallelComp only opens a pool
% and the fits run one after another either way
if fitS.useParallelComp
    if sum(size(gcp("nocreate"))) == 0
        parpool('Processes', 4)
    end
else
    parallelSettings                 = parallel.Settings;
    poolAutoCreateBefore             = parallelSettings.Pool.AutoCreate;
    parallelSettings.Pool.AutoCreate = false;
end

%% Initialise the output structure
% Initialise the data and summary fields (one column or cell per iteration)
outS.OneRunDatas = OneRunDatas;
outS.agediffsWbC = NaN(length(dtBinEdges_yr)-1, numruns);
outS.weightedC   = cell(1, numruns);
outS.numCpairs   = NaN(1, numruns);
outS.sedLength   = NaN(1, numruns);
outS.sedTimeSpan = NaN(1, numruns);

% Record the per-method choice with the results, since fitS is saved once
% for all methods
outS.Lin2014AgeFiltering = applyLin2014AgeFilter;

% Initialise each distribution's pdf on the grid x (.nSR) and on log(x)
% (.lnSR), with their mean and variance, and its fit results (.fits)

% --- MLN
outS.MLN.nSR.x    = x;
outS.MLN.nSR.px   = NaN(numX, numruns);
outS.MLN.nSR.mu   = NaN(1, numruns);
outS.MLN.nSR.var  = NaN(1, numruns);
outS.MLN.lnSR.x   = logx;
outS.MLN.lnSR.px  = NaN(numX, numruns);
outS.MLN.lnSR.mu  = NaN(1, numruns);
outS.MLN.lnSR.var = NaN(1, numruns);
outS.MLN.fits     = cell(1, numruns);

% --- invGam
outS.invGam.nSR.x    = x;
outS.invGam.nSR.px   = NaN(numX, numruns);
outS.invGam.nSR.mu   = NaN(1, numruns);
outS.invGam.nSR.var  = NaN(1, numruns);
outS.invGam.lnSR.x   = logx;
outS.invGam.lnSR.px  = NaN(numX, numruns);
outS.invGam.lnSR.mu  = NaN(1, numruns);
outS.invGam.lnSR.var = NaN(1, numruns);
outS.invGam.fits     = cell(1, numruns);

% --- Gam
outS.Gam.nSR.x    = x;
outS.Gam.nSR.px   = NaN(numX, numruns);
outS.Gam.nSR.mu   = NaN(1, numruns);
outS.Gam.nSR.var  = NaN(1, numruns);
outS.Gam.lnSR.x   = logx;
outS.Gam.lnSR.px  = NaN(numX, numruns);
outS.Gam.lnSR.mu  = NaN(1, numruns);
outS.Gam.lnSR.var = NaN(1, numruns);
outS.Gam.fits     = cell(1, numruns);

% --- LN
outS.LN.nSR.x    = x;
outS.LN.nSR.px   = NaN(numX, numruns);
outS.LN.nSR.mu   = NaN(1, numruns);
outS.LN.nSR.var  = NaN(1, numruns);
outS.LN.lnSR.x   = logx;
outS.LN.lnSR.px  = NaN(numX, numruns);
outS.LN.lnSR.mu  = NaN(1, numruns);
outS.LN.lnSR.var = NaN(1, numruns);
outS.LN.fits     = cell(1, numruns);

% Initialise the chi-squared results (h and p have one row per pdf and one
% column per iteration)
h        = NaN(4, numruns);
p        = NaN(4, numruns);
chiStats = cell(numruns, 1);

%% Fit the distributions to each iteration's data
for i = 1:numruns

    % --- Prepare the data
    % Get this iteration's data
    OneRunData = OneRunDatas{i};

    % Remove negative SRs, keeping the NaN run-start columns (ideally there
    % shouldn't be any, but BMode can produce them)
    if sum(OneRunData(1,:) < 0) ~= 0
        warning("There are negative sed rates!... being removed")
        zeroOrPositiveLog = OneRunData(1,:) >= 0 | isnan(OneRunData(1,:));
        OneRunData        = OneRunData(:, zeroOrPositiveLog);
    end

    % Remove SRs of 0 (ideally there shouldn't be any)
    if sum(OneRunData(1,:) == 0) ~= 0
        warning("There are SRs of 0!... being removed.")
        notZeroLog = OneRunData(1,:) ~= 0;
        OneRunData = OneRunData(:, notZeroLog);
    end

    % If wanted, use absolute SR (cm/kyr) instead of nSR. Row 2 is dx in cm
    % and row 3 is dt in yr, so the ratio is multiplied by 1000 to give
    % cm/kyr
    if fitS.non_normalized_SR
        dx_cm = OneRunData(2,:);
        dt_yr = OneRunData(3,:);
        SRs   = dx_cm./dt_yr.*1000;

        % Keep any NaN nSRs as NaN, and put the SRs in row 1
        nanLog          = isnan(OneRunData(1,:));
        SRs(nanLog)     = NaN;
        OneRunData(1,:) = SRs;
    end

    % Apply the age filter as done by Lin2014 if desired: keep the pairs
    % whose dt_yr is at least the lower value of fitS.Lin2014AgeFilter and
    % below the upper value (the lower bound is inclusive, as in ARfitdists,
    % so that pairs with a dt of exactly 500 yr, which whole-year ages give,
    % are kept), and any NaN nSRs
    if applyLin2014AgeFilter
        maxDt_yr       = max(fitS.Lin2014AgeFilter);
        minDt_yr       = min(fitS.Lin2014AgeFilter);
        inAgeFilterLog = OneRunData(3,:) < maxDt_yr & OneRunData(3,:) >= minDt_yr;
        nanLog         = isnan(OneRunData(1,:));
        lin2014KeepLog = inAgeFilterLog | nanLog;
        OneRunData     = OneRunData(:, lin2014KeepLog);
    end

    % Set up nSR counts and weighting vectors without any NaN nSRs (nSRs are
    % SRs in cm/kyr if fitS.non_normalized_SR)
    notNanLog = ~isnan(OneRunData(1,:));
    nSRs      = OneRunData(1, notNanLog);
    dx_cm     = OneRunData(2, notNanLog);
    dt_yr     = OneRunData(3, notNanLog);

    % Get the replication settings for the individual-run fits
    weightRepDP       = fitS.OneRun.weightRepDP;
    weightRepInflator = fitS.OneRun.weightRepInflator;

    % Set up the weights according to the weighting choice, and weight by
    % replicating the data according to the weights (not for "none")
    if fitS.weighting == "none"
        weights   = ones(size(nSRs));
        weightedC = nSRs;
    elseif fitS.weighting == "depth"
        weights   = dx_cm;   % Weight by depth difference
        weightedC = makeWeightedReplicates(nSRs, weights, weightRepDP, weightRepInflator);
    elseif fitS.weighting == "age"
        weights   = dt_yr;   % Weight by time span
        weightedC = makeWeightedReplicates(nSRs, weights, weightRepDP, weightRepInflator);
    else
        error("Weighting type not recognized, must be 'none', 'depth', or 'age'")
    end

    % Calculate the summed dx (cm) and dt (yr) of this iteration's data (over
    % the pairs, not the run-start columns), the number of nSR values, and
    % the weighted bin counts of dt_yr
    sedLength        = sum(dx_cm);
    sedTimeSpan      = sum(dt_yr);
    numSRcalcs       = numel(nSRs);
    dtWeightedCounts = makeWeightedBinCounts(dt_yr, weights, dtBinEdges_yr);

    % Resample data back to the original number of nSR values if desired
    % (removes the inflation from the weighting)
    if fitS.resampleData
        rng(1);
        weightedC = randsample(weightedC, numSRcalcs, false);
    end

    % --- Fit the distributions
    % Fit the four distributions, trying again if a fit fails, or store NaN
    % placeholders if fitting is switched off (fitS.fitDists)
    skipIteration = false; % This will change to true if 10 iterations fail, exiting while loop

    if fitS.fitDists
        % Make sure the data are a row vector for the fitting functions
        if size(weightedC, 1) > size(weightedC, 2)
            weightedC = weightedC';
        end

        % Fit the 2-component MLN, the invGam, the Gam and the LN (a
        % 1-component MLN), each giving its pdf on the grid x
        fitsSucceeded = false;
        numFailedFits = 0;
        while ~fitsSucceeded && ~skipIteration
            try
                [mixLogNormPDF, ~, MLNfitInfo] = fitMixLogNorm(weightedC, x, 2, ...
                    fitS.OneRun.MLNReps, numSRcalcs);
                [invGamPDF, ~, invGamFitInfo]  = fitInvGamma(weightedC, x, numSRcalcs);
                [GamPDF, GamFitInfo]           = fitGamma(weightedC, x, numSRcalcs);
                [logNormPDF, ~, LNfitInfo]     = fitMixLogNorm(weightedC, x, 1, ...
                    fitS.OneRun.MLNReps, numSRcalcs);

                % Replace any NaNs in the pdfs with 0 (think this is a bug in
                % fitters?)
                mixLogNormPDF(isnan(mixLogNormPDF)) = 0;
                invGamPDF(isnan(invGamPDF))         = 0;
                GamPDF(isnan(GamPDF))               = 0;
                logNormPDF(isnan(logNormPDF))       = 0;

                fitsSucceeded = true;
            catch
                % fitMixLogNorm can fail with an ill-conditioned covariance;
                % retrying with new random initial conditions often succeeds.
                % If it fails more than 10 times, skip this iteration
                numFailedFits = numFailedFits + 1;
                disp("number of failed fittings = " + num2str(numFailedFits))
                if numFailedFits > 10
                    skipIteration = true;
                    [MLNfitInfo, invGamFitInfo, GamFitInfo, LNfitInfo, ...
                        mixLogNormPDF, invGamPDF, GamPDF, logNormPDF] = makeDummyFitResult(numX);
                end
            end

            % Report progress every 50 iterations
            if mod(i, 50) == 0
                disp("Fit " + num2str(i) + " samplings")
            end
        end
    else
        % Distribution fitting is switched off, so store NaN placeholders
        [MLNfitInfo, invGamFitInfo, GamFitInfo, LNfitInfo, ...
            mixLogNormPDF, invGamPDF, GamPDF, logNormPDF] = makeDummyFitResult(numX);
    end

    % Report a skipped iteration
    if skipIteration
        disp("Skipped an iteration in IRfitdists")
    end

    % --- Build the pdf structs
    % Make each pdf into a struct on the grid x and one on log(x), with its
    % mean and variance (NaN for the NaN pdfs of a skipped or switched-off
    % fit)
    MLNVEC       = buildVEC(x, mixLogNormPDF(:,2), 5, "2 Component Mixed Log Normal");
    logMLNVEC    = buildLogVEC(MLNVEC);
    invGamVEC    = buildVEC(x, invGamPDF,          2, "Inverse Gamma");
    logInvGamVEC = buildLogVEC(invGamVEC);
    GamVEC       = buildVEC(x, GamPDF,             2, "Gamma");
    logGamVEC    = buildLogVEC(GamVEC);
    LNVEC        = buildVEC(x, logNormPDF(:,2),    2, "Log Normal");
    logLNVEC     = buildLogVEC(LNVEC);

    % --- Run the chi-squared goodness-of-fit tests if desired
    % Test log(weightedC) against each log pdf, starting from 20 bins and
    % plotting the tests of the first three iterations; otherwise store NaN
    % placeholders
    % TODO: fitS.dispChi2 is set here, so fitData's setting has no effect on
    % these tests, and the 20 initial bins could come from fitS
    % (fitS.chi2binN is not used here)
    if ~skipIteration && fitS.run_chi2gof
        pdfs          = {logLNVEC; logMLNVEC; logGamVEC; logInvGamVEC};
        fitS.dispChi2 = (i < 4);
        [run_h, run_p, run_chiStats] = chi2_dataVStwopdfVECs(log(weightedC), numSRcalcs, 20, ...
            pdfs, fitS);
    else
        [run_h, run_p, run_chiStats] = makeDummyChi2(numPdfs);
    end

    % --- Store the results of this iteration
    outS.agediffsWbC(:,i) = dtWeightedCounts';
    outS.weightedC{i}     = weightedC;
    outS.numCpairs(i)     = numSRcalcs;
    outS.sedLength(i)     = sedLength;
    outS.sedTimeSpan(i)   = sedTimeSpan;

    outS.MLN.fits{i}      = MLNfitInfo;
    outS.MLN.nSR.px(:,i)  = MLNVEC.px;
    outS.MLN.nSR.mu(i)    = MLNVEC.mu;
    outS.MLN.nSR.var(i)   = MLNVEC.var;
    outS.MLN.lnSR.px(:,i) = logMLNVEC.px;
    outS.MLN.lnSR.mu(i)   = logMLNVEC.mu;
    outS.MLN.lnSR.var(i)  = logMLNVEC.var;

    outS.invGam.fits{i}      = invGamFitInfo;
    outS.invGam.nSR.px(:,i)  = invGamVEC.px;
    outS.invGam.nSR.mu(i)    = invGamVEC.mu;
    outS.invGam.nSR.var(i)   = invGamVEC.var;
    outS.invGam.lnSR.px(:,i) = logInvGamVEC.px;
    outS.invGam.lnSR.mu(i)   = logInvGamVEC.mu;
    outS.invGam.lnSR.var(i)  = logInvGamVEC.var;

    outS.Gam.fits{i}      = GamFitInfo;
    outS.Gam.nSR.px(:,i)  = GamVEC.px;
    outS.Gam.nSR.mu(i)    = GamVEC.mu;
    outS.Gam.nSR.var(i)   = GamVEC.var;
    outS.Gam.lnSR.px(:,i) = logGamVEC.px;
    outS.Gam.lnSR.mu(i)   = logGamVEC.mu;
    outS.Gam.lnSR.var(i)  = logGamVEC.var;

    outS.LN.fits{i}      = LNfitInfo;
    outS.LN.nSR.px(:,i)  = LNVEC.px;
    outS.LN.nSR.mu(i)    = LNVEC.mu;
    outS.LN.nSR.var(i)   = LNVEC.var;
    outS.LN.lnSR.px(:,i) = logLNVEC.px;
    outS.LN.lnSR.mu(i)   = logLNVEC.mu;
    outS.LN.lnSR.var(i)  = logLNVEC.var;

    h(:,i)      = run_h';
    p(:,i)      = run_p';
    chiStats{i} = run_chiStats;

end

%% Restore the parallel pool setting
% Let MATLAB open a pool automatically again, if it could before
if ~fitS.useParallelComp
    parallelSettings.Pool.AutoCreate = poolAutoCreateBefore;
end

%% Put the chi-squared results into one table per distribution
% The pdfs were tested in the order LN, MLN, Gam, invGam (the pdfs cell in
% the fitting loop), so that is the row each distribution's results are in
outS.MLN.chiStats    = buildChi2Table(chiStats, h(2,:)', p(2,:)', 2);
outS.invGam.chiStats = buildChi2Table(chiStats, h(4,:)', p(4,:)', 4);
outS.Gam.chiStats    = buildChi2Table(chiStats, h(3,:)', p(3,:)', 3);
outS.LN.chiStats     = buildChi2Table(chiStats, h(1,:)', p(1,:)', 1);

end

%% Local functions
% Helpers for the fitting loop above: NaN placeholders for a skipped fit or
% test, the pdf structs, and the chi-squared tables

function[MLNfitInfo, invGamFitInfo, GamFitInfo, LNfitInfo, ...
    mixLogNormPDF, invGamPDF, GamPDF, logNormPDF] = makeDummyFitResult(numX)
% Returns NaN-filled fit results and pdfs, for when fitting is switched off
% (fitS.fitDists = false) or has failed more than 10 times
MLNfitInfo    = struct('NumParams', NaN, 'mu', NaN, 'Sigma', NaN, ...
                       'ComponentProportion', NaN, 'NegativeLogLikelihood', NaN, ...
                       'BIC', NaN, 'BICtaeheefix', NaN);
invGamFitInfo = struct('NumParams', NaN, 'alpha', NaN, 'beta', NaN, ...
                       'NegativeLogLikelihood', NaN, 'BIC', NaN, 'BICtaeheefix', NaN);
GamFitInfo    = struct('NumParams', NaN, 'alpha', NaN, 'beta', NaN, ...
                       'NegativeLogLikelihood', NaN, 'BIC', NaN, 'BICtaeheefix', NaN);
LNfitInfo     = struct('NumParams', NaN, 'mu', NaN, 'Sigma', NaN, ...
                       'NegativeLogLikelihood', NaN, 'BIC', NaN, 'BICtaeheefix', NaN);
mixLogNormPDF = NaN(numX, 2);
invGamPDF     = NaN(numX, 1);
GamPDF        = NaN(numX, 1);
logNormPDF    = NaN(numX, 2);
end

function[pdfVEC] = buildVEC(x, px, numParams, pdfName)
% Makes a pdf struct and calculates its mean and variance (NaN if px is NaN,
% from a skipped fit)
pdfVEC.x         = x';
pdfVEC.px        = px;
pdfVEC.numParams = numParams;
pdfVEC.pdfName   = pdfName;
[pdfVEC.mu, pdfVEC.var] = muVarPDFVec(pdfVEC);
end

function[logPdfVEC] = buildLogVEC(pdfVEC)
% Transforms a pdf struct from nSR to log(nSR) and calculates its mean and
% variance there
[logPdfVEC.x, logPdfVEC.px] = px_to_pfx(pdfVEC.x, pdfVEC.px, @log);
logPdfVEC.numParams         = pdfVEC.numParams;
logPdfVEC.pdfName           = pdfVEC.pdfName;
[logPdfVEC.mu, logPdfVEC.var] = muVarPDFVec(logPdfVEC);
end

function[run_h, run_p, run_chiStats] = makeDummyChi2(numPdfs)
% Returns NaN-filled chi-squared results for one iteration, for when the
% tests are skipped or switched off
run_h        = nan(1, numPdfs);
run_p        = nan(1, numPdfs);
run_chiStats = cell(1, numPdfs);
for i = 1:numPdfs
    run_chiStats{i} = struct('chi2stat', NaN, 'df', NaN, 'edges', NaN, 'O', NaN, 'E', NaN);
end
end

function[chiStatsT] = buildChi2Table(chiStats, h, p, distIdx)
% Puts one distribution's chi-squared results from every iteration into a
% table (one row per iteration). distIdx is the position of the
% distribution's results in each iteration's cell, the same as its position
% in the pdfs cell in the fitting loop
numRuns  = length(chiStats);
chi2stat = cell(numRuns, 1);
df       = NaN(numRuns, 1);
edges    = cell(numRuns, 1);
O        = cell(numRuns, 1);
E        = cell(numRuns, 1);
for i = 1:numRuns
    chi2stat{i} = chiStats{i}{distIdx}.chi2stat;
    df(i)       = chiStats{i}{distIdx}.df;
    edges{i}    = chiStats{i}{distIdx}.edges;
    O{i}        = chiStats{i}{distIdx}.O;
    E{i}        = chiStats{i}{distIdx}.E;
end
chiStatsT = table(chi2stat, h, p, df, edges, O, E, ...
    'VariableNames', ["chi2stat","h","p","df","edges","O","E"]);
end
