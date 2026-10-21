%% BICeff_nSensitivity - BICeff of the four families against the effective number of pairs (Reviewer 1)
%
% Reviewer 1 objected that BICeff corrects the sample-size term for the
% replication and pooling, but that autocorrelation within cores and
% dependence among the weighted replicates bias the likelihood, so the
% ranking of the families could be wrong. The log-likelihood of the
% replicated vector the fitters saw (weightedC, numWP points) is a weighted
% log-likelihood of the n_RP = numCpairs radiocarbon pairs, and BICeff
% rescales it to the per-point average times n_RP, so the only open
% question is which effective n to use. The fitted parameters don't change
% with n, so BICeff(n) = -2*n*lbar + k*ln(n), with lbar = lnL/numWP the
% per-point average log-likelihood of the fit: the ranking among the k = 2
% families (LN, Gam, invGam) doesn't depend on n at all, and LN vs MLN
% (k = 5) changes sign at a crossover n*. This script gets lbar from the
% manuscript fit file, checks that BICeff(n_RP) reproduces the stored
% BICeff, and gives BICeff and the ranking at the benchmarks: the number of
% fitted cores, the Kish effective n of the depth weights of one
% realisation (autocorrelation and unequal weights both reduce the
% effective n), multiples of n_RP, and numWP; then n* itself. The per-run
% disclosure (how often MLN beats LN in a single-run fit) comes from the
% individual-run fits (<method>IR) in the same file.
%
% The script cds to the repository root because d.dataFile is a path
% relative to it.

cd(fileparts(fileparts(mfilename('fullpath'))))
addpath('Functions')

%% Settings
fitFile         = "Results/results_QSR/dataT_QSR_fit_depthweight.mat"; % Manuscript fit
runTag          = "QSR";               % Tag of the calcData run, used in the saved file names
methodNames     = ["BMedianLin", "BSampLin", "RSR0", "RSR500", "RSR1000"]; % The manuscript's
                                       %   Table 2 methods; any other fit in fitFile can be added
distNames       = ["LN", "MLN", "Gam", "invGam"]; % Field names in the fit structs
numParams       = [2 5 2 2];           % Number of fitted parameters (k) of each family in distNames
nMultiples      = [0.5 1 2 5 10];      % Multiples of n_RP (numCpairs) tabulated as benchmarks
nSweep          = logspace(1, 8, 500); % n values for the figure
nCrossoverGrid  = logspace(0, 8, 1e6); % Fine grid of n searched for the LN-MLN crossover (spacing
                                       %   about 0.002 % of n)
BICeffTolerance = 1e-6;                % Largest |recalculated - stored BICeff| allowed (rounding only)
figYLim         = [-40 40];            % y-axis limits of Fig 1 (BIC units); BIC_LN - BIC_MLN grows
                                       %   linearly with n, so without limits the crossover is flattened
resultsFolder   = "Results/review1_nSensitivity"; % Where the tables are saved
useKish         = true;                % If true, loads the calcData file (d.dataFile, ~19 s) for the
                                       %   Kish effective n of the depth weights

%% Load the fit results and keep what is needed of each method
% d is about 1.7 GB, so the numbers needed are copied into plain arrays and
% d is cleared afterwards
disp("Loading " + fitFile)
load(fitFile, "d")

% Initialise the arrays (one row per method, one column per family)
numMethods      = numel(methodNames);
numDists        = numel(distNames);
NLL             = NaN(numMethods, numDists);  % NLL of all the weighted points the fitters saw
storedBICeff    = NaN(numMethods, numDists);  % BICeff as saved by the fitters (fitInfo.BICtaeheefix)
storedNumParams = NaN(numMethods, numDists);  % k as saved by the fitters (fitInfo.NumParams)
numWP           = NaN(numMethods, 1);         % Number of weighted points fitted (numel(weightedC))
numCpairs       = NaN(numMethods, 1);         % n_RP: number of nSR pairs per run
lin2014Log      = false(numMethods, 1);       % Whether the fit had the Lin2014 age filter
hasMethodLog    = false(numMethods, 1);
hasIRLog        = false(numMethods, 1);
IRLNminusMLN    = cell(numMethods, 1);        % BICeff_LN - BICeff_MLN of each individual-run fit

for i = 1:numMethods
    % Skip a method that isn't in the fit file
    if ~isfield(d.S1, methodNames(i))
        warning(methodNames(i) + " is not in " + fitFile + ", so it is skipped")
        continue
    end
    hasMethodLog(i) = true;
    methodFit       = d.S1.(methodNames(i));

    % Get the number of pairs per run, the number of points fitted and the
    % age-filter choice
    numCpairs(i)  = methodFit.numCpairs;
    numWP(i)      = numel(methodFit.weightedC);
    lin2014Log(i) = methodFit.Lin2014AgeFiltering;

    % Get the NLL, the stored BICeff and k of each family
    for j = 1:numDists
        distFitInfo          = methodFit.(distNames(j)).fitInfo;
        NLL(i,j)             = distFitInfo.NegativeLogLikelihood;
        storedBICeff(i,j)    = distFitInfo.BICtaeheefix;
        storedNumParams(i,j) = distFitInfo.NumParams;
    end

    % Find the individual-run fits of this method, if there are any
    IRName = methodNames(i) + "IR";
    if ~isfield(d.S1, IRName)
        disp(methodNames(i) + " has no individual-run fits (" + IRName + ")")
        continue
    end
    hasIRLog(i) = true;

    % Get BICeff_LN - BICeff_MLN of each run (each run's BICeff uses that
    % run's own number of pairs, IR.numCpairs)
    numIRruns      = numel(d.S1.(IRName).LN.fits);
    run_LNminusMLN = NaN(numIRruns, 1);
    for iRun = 1:numIRruns
        run_BICeffLN         = d.S1.(IRName).LN.fits{iRun}.BICtaeheefix;
        run_BICeffMLN        = d.S1.(IRName).MLN.fits{iRun}.BICtaeheefix;
        run_LNminusMLN(iRun) = run_BICeffLN - run_BICeffMLN;
    end
    IRLNminusMLN{i} = run_LNminusMLN;
end

% Get the fitted cores, the calcData file and the fit settings
numCores        = d.S1.numCores;
chosenCoreNames = d.S1.chosenCoreNames;
dataFile        = d.dataFile;
fitS            = d.S1.fitS;
clear d methodFit

% Check that k in the settings is the k the fitters used
for j = 1:numDists
    kMismatchLog = hasMethodLog & storedNumParams(:,j) ~= numParams(j);
    if any(kMismatchLog)
        error("numParams for " + distNames(j) + " doesn't match fitInfo.NumParams in the fit file")
    end
end

%% Calculate the per-point average log-likelihood and check it reproduces the stored BICeff
% lbar = lnL/numWP = -NLL/numWP. The fitters divide NLL by the replication
% factor numWP/numCpairs and add k*ln(numCpairs), so BICeff(n) = -2*n*lbar +
% k*ln(n) at n = numCpairs must give back fitInfo.BICtaeheefix
meanLogLik = -NLL./numWP;

% Recalculate BICeff at n = numCpairs
recalcBICeff = NaN(numMethods, numDists);
for j = 1:numDists
    recalcBICeff(:,j) = -2.*numCpairs.*meanLogLik(:,j) + numParams(j).*log(numCpairs);
end

% Stop if it doesn't match: every number below rests on this
BICeffDiffs   = abs(recalcBICeff - storedBICeff);
maxBICeffDiff = max(BICeffDiffs(:), [], 'omitnan');
disp("Max |BICeff(n_RP) recalculated from lbar - fitInfo.BICtaeheefix| = " + maxBICeffDiff)
if maxBICeffDiff > BICeffTolerance
    error("BICeff at n = numCpairs doesn't reproduce fitInfo.BICtaeheefix: either " + ...
        "fitInfo.NegativeLogLikelihood is not the NLL of all of weightedC, or numCpairs " + ...
        "is not the n the fitters used")
end

%% Calculate the Kish effective n of the depth weights of one realisation
% Kish n = (sum w)^2 / sum(w^2), with w the dx_cm of each pair (row 2 of
% the pair columns of the nSR matrix), pooled over the fitted cores. One
% realisation only: run 1 of the Monte Carlo methods, the single
% realisation of BMedian. The pairs kept are the ones ARfitdists fits:
% positive nSR, and for the <method>Lin fits dt_yr strictly between the
% two values of fitS.Lin2014AgeFilter. Not sure whether run 1 is typical;
% kishNumPairs is printed next to n_RP (the mean over runs) as a check
kishN        = NaN(numMethods, 1);
kishNumPairs = NaN(numMethods, 1);   % Number of pairs in the realisation used for Kish n

if useKish
    % Warn where the weights here wouldn't be the ones fitted
    if fitS.weighting ~= "depth"
        warning("The fits were weighted by " + fitS.weighting + ", but Kish n uses the depth weights")
    end

    % Load the calcData results
    disp("Loading " + dataFile + " for the Kish effective n")
    load(dataFile, "dataT")

    % Find the fitted cores
    fittedCoreLog = ismember(string(dataT.cores), chosenCoreNames);
    if sum(fittedCoreLog) ~= numCores
        warning("Found " + sum(fittedCoreLog) + " of the " + numCores + " fitted cores in " + dataFile)
    end

    % Set up the Lin2014 age filter limits (yr)
    maxDt_yr = max(fitS.Lin2014AgeFilter);
    minDt_yr = min(fitS.Lin2014AgeFilter);

    for i = 1:numMethods
        if ~hasMethodLog(i)
            continue
        end

        % Find the dataT column: the <method>Lin fits are fits of <method>
        dataTColName = methodNames(i);
        if lin2014Log(i)
            dataTColName = extractBefore(methodNames(i), "Lin");
        end
        if ~ismember(dataTColName, dataT.Properties.VariableNames)
            warning(dataTColName + " is not in dataT, so " + methodNames(i) + " has no Kish n")
            continue
        end

        % Pool the depth weights of run 1 of each fitted core
        pooledWeights = [];
        for j = 1:height(dataT)
            core_nSRmatrix = dataT.(dataTColName){j};
            if ~fittedCoreLog(j) || isempty(core_nSRmatrix)
                continue
            end

            % Get run 1 (header column [NaN; depth; age], then [nSR; dx_cm; dt_yr] per pair)
            run_nSRmatrix = selectRuns(core_nSRmatrix, 1);

            % Keep the pairs with positive nSR (this also drops the NaN header column)
            pairLog = run_nSRmatrix(1,:) > 0;

            % Apply the Lin2014 age filter if the fit had it
            if lin2014Log(i)
                inAgeFilterLog = run_nSRmatrix(3,:) < maxDt_yr & run_nSRmatrix(3,:) > minDt_yr;
                pairLog        = pairLog & inAgeFilterLog;
            end

            % Add this core's weights
            run_dx_cm     = run_nSRmatrix(2, pairLog);
            pooledWeights = [pooledWeights, run_dx_cm]; %#ok<AGROW>
        end

        % Calculate the Kish effective n
        kishN(i)        = sum(pooledWeights).^2./sum(pooledWeights.^2);
        kishNumPairs(i) = numel(pooledWeights);
        disp(methodNames(i) + ": n_RP = " + compose("%.0f", numCpairs(i)) + ", pairs in run 1 = " + ...
            kishNumPairs(i) + ", Kish n = " + compose("%.0f", kishN(i)))
    end
    clear dataT
end

%% Calculate BICeff of every family at each benchmark n
% Benchmarks: the number of fitted cores, the Kish n, n_RP times each of
% nMultiples (the x1 row is n_RP itself), and numWP (the replicated counts)
benchmarkNames = ["numCores", "kishN", "nRP x" + nMultiples, "nWP"];
numBenchmarks  = numel(benchmarkNames);
numRows        = numMethods.*numBenchmarks;

% Initialise the table columns
methodCol    = strings(numRows, 1);
benchmarkCol = strings(numRows, 1);
nCol         = NaN(numRows, 1);
BICeffCols   = NaN(numRows, numDists);
rankingCol   = strings(numRows, 1);
LNminusMLN   = NaN(numRows, 1);

% Find the LN and MLN columns
LNIdx  = find(distNames == "LN");
MLNIdx = find(distNames == "MLN");

iRow = 0;
for i = 1:numMethods
    if ~hasMethodLog(i)
        continue
    end

    % Set up this method's benchmark n values, in the order of benchmarkNames
    benchmarkN = [numCores, kishN(i), nMultiples.*numCpairs(i), numWP(i)];

    for iBench = 1:numBenchmarks
        % Record which method, benchmark and n this row is
        iRow = iRow + 1;
        benchN = benchmarkN(iBench);
        methodCol(iRow)    = methodNames(i);
        benchmarkCol(iRow) = benchmarkNames(iBench);
        nCol(iRow)         = benchN;

        % Calculate BICeff of each family at this n
        bench_BICeff = -2.*benchN.*meanLogLik(i,:) + numParams.*log(benchN);
        BICeffCols(iRow,:) = bench_BICeff;
        LNminusMLN(iRow)   = bench_BICeff(LNIdx) - bench_BICeff(MLNIdx);

        % Rank the families (lowest BICeff first); no ranking without an n
        if isnan(benchN)
            continue
        end
        [~, distOrder]   = sort(bench_BICeff);
        rankingCol(iRow) = strjoin(distNames(distOrder), " > ");
    end
end

% Put the columns into a table, leaving out the rows of skipped methods
filledLog  = methodCol ~= "";
benchmarkT = table(methodCol, benchmarkCol, nCol, 'VariableNames', ["method", "benchmarkName", "n"]);
for j = 1:numDists
    benchmarkT.("BICeff_" + distNames(j)) = BICeffCols(:,j);
end
benchmarkT.ranking    = rankingCol;
benchmarkT.LNminusMLN = LNminusMLN;
benchmarkT            = benchmarkT(filledLog,:);

%% Find the n at which MLN overtakes LN
% BIC_LN - BIC_MLN = 2*n*(lbar_MLN - lbar_LN) - (k_MLN - k_LN)*ln(n), and LN
% is preferred where it is negative. It is convex in n, so LN wins on one
% interval at most: for gainMLN = lbar_MLN - lbar_LN > 0 the difference
% falls from 2*gainMLN at n = 1 to a minimum at n = (k_MLN -
% k_LN)/(2*gainMLN) and rises again, so it has two roots, a small one and
% the crossover n* above which MLN wins. I search nCrossoverGrid for the
% interval rather than solving the equation, so the ends are good to the
% grid spacing. lowN is the smallest n at which LN wins, crossoverN (n*)
% the first n above the interval at which MLN wins
extraParams = numParams(MLNIdx) - numParams(LNIdx);
gainMLN     = meanLogLik(:,MLNIdx) - meanLogLik(:,LNIdx);   % Per-point log-likelihood MLN gains
lowN        = NaN(numMethods, 1);
crossoverN  = NaN(numMethods, 1);

for i = 1:numMethods
    if ~hasMethodLog(i)
        continue
    end

    % Calculate BIC_LN - BIC_MLN on the grid and find where LN wins
    grid_LNminusMLN = 2.*nCrossoverGrid.*gainMLN(i) - extraParams.*log(nCrossoverGrid);
    LNWinsLog       = grid_LNminusMLN < 0;
    LNWinsIdx       = find(LNWinsLog);

    % MLN wins at every n on the grid: leave both NaN
    if isempty(LNWinsIdx)
        disp(methodNames(i) + ": MLN is preferred at every n on nCrossoverGrid")
        continue
    end
    lowN(i) = nCrossoverGrid(LNWinsIdx(1));

    % LN still wins at the top of the grid: no crossover within it
    if LNWinsIdx(end) == numel(nCrossoverGrid)
        crossoverN(i) = Inf;
        continue
    end
    crossoverN(i) = nCrossoverGrid(LNWinsIdx(end) + 1);
end

%% Summarise the per-run disclosure from the individual-run fits
% For each run: BICeff_LN - BICeff_MLN; positive means MLN is preferred in
% that run. Fraction of runs where MLN is preferred, and the median and the
% 5-95 % range of the difference (runs with a failed fit, NaN, left out)
fracRunsMLN   = NaN(numMethods, 1);
IRmedian      = NaN(numMethods, 1);
IRpct5        = NaN(numMethods, 1);
IRpct95       = NaN(numMethods, 1);
numIRrunsUsed = NaN(numMethods, 1);

for i = 1:numMethods
    if ~hasIRLog(i)
        continue
    end

    % Get the runs with a finite difference
    run_LNminusMLN = IRLNminusMLN{i};
    finiteRunLog   = isfinite(run_LNminusMLN);
    run_LNminusMLN = run_LNminusMLN(finiteRunLog);

    % Calculate the fraction where MLN is preferred, and the spread
    numIRrunsUsed(i) = numel(run_LNminusMLN);
    fracRunsMLN(i)   = sum(run_LNminusMLN > 0)./numIRrunsUsed(i);
    IRmedian(i)      = median(run_LNminusMLN);
    IRpct5(i)        = prctile(run_LNminusMLN, 5);
    IRpct95(i)       = prctile(run_LNminusMLN, 95);
end

%% Rank the families with the same k as LN (independent of n)
% With equal k the k*ln(n) terms cancel, so these families are ranked by
% lbar alone (highest first) at every n
sameKLog     = numParams == numParams(LNIdx);
sameKNames   = distNames(sameKLog);
sameKRanking = strings(numMethods, 1);
for i = 1:numMethods
    if ~hasMethodLog(i)
        continue
    end
    [~, sameKOrder] = sort(meanLogLik(i, sameKLog), 'descend');
    sameKRanking(i) = strjoin(sameKNames(sameKOrder), " > ");
end

%% Save the tables
% One row per method x benchmark, and one row per method for the crossover
% and the per-run numbers
if ~isfolder(resultsFolder)
    mkdir(resultsFolder)
end
writetable(benchmarkT, fullfile(resultsFolder, "BICeff_nSensitivity_" + runTag + ".csv"))

% Set up the crossover and per-run table
crossoverT = table(methodNames', numCores.*ones(numMethods, 1), kishN, kishNumPairs, numCpairs, ...
    numWP, gainMLN, lowN, crossoverN, crossoverN./numCpairs, sameKRanking, ...
    'VariableNames', ["method", "numCores", "kishN", "kishNumPairs", "nRP", "nWP", "gainMLN", ...
    "lowN", "crossoverN", "crossoverOverNRP", "sameKRanking"]);
crossoverT.numIRruns          = numIRrunsUsed;
crossoverT.IRfracMLNwins      = fracRunsMLN;
crossoverT.IRmedianLNminusMLN = IRmedian;
crossoverT.IRpct5LNminusMLN   = IRpct5;
crossoverT.IRpct95LNminusMLN  = IRpct95;
crossoverT                    = crossoverT(hasMethodLog,:);
writetable(crossoverT, fullfile(resultsFolder, "BICeff_crossover_" + runTag + ".csv"))
disp("Saved the tables in " + resultsFolder)

%% Print a summary
% One line per method: n_RP, Kish n, n*, the interval where LN wins, which
% family wins at n_RP, and the k = 2 ranking; then the per-run numbers
disp(" ")
disp("BICeff against n (BIC_LN - BIC_MLN = 2*n*gainMLN - " + extraParams + "*ln(n)):")
for i = 1:numMethods
    if ~hasMethodLog(i)
        continue
    end

    % Find which of LN and MLN wins at n_RP
    nRP_LNminusMLN = 2.*numCpairs(i).*gainMLN(i) - extraParams.*log(numCpairs(i));
    winnerAtNRP    = "LN";
    if nRP_LNminusMLN > 0
        winnerAtNRP = "MLN";
    end

    disp(methodNames(i) + ": n_RP = " + compose("%.0f", numCpairs(i)) + ", Kish n = " + ...
        compose("%.0f", kishN(i)) + ", n* = " + compose("%.0f", crossoverN(i)) + " (" + ...
        compose("%.2f", crossoverN(i)./numCpairs(i)) + " x n_RP), LN wins for " + ...
        compose("%.1f", lowN(i)) + " < n < n*, " + winnerAtNRP + " wins at n_RP; " + ...
        "k=" + numParams(LNIdx) + " ranking (n-independent): " + sameKRanking(i))
end

disp(" ")
disp("Individual-run fits (BICeff_LN - BICeff_MLN per run; > 0 means MLN preferred):")
for i = 1:numMethods
    if ~hasIRLog(i)
        disp(methodNames(i) + ": no individual-run fits")
        continue
    end
    disp(methodNames(i) + "IR: MLN preferred in " + compose("%.3f", fracRunsMLN(i)) + " of " + ...
        numIRrunsUsed(i) + " runs, median " + compose("%.1f", IRmedian(i)) + ", 5-95 % [" + ...
        compose("%.1f", IRpct5(i)) + ", " + compose("%.1f", IRpct95(i)) + "]")
end

%% Fig 1: BIC_LN - BIC_MLN against n, one panel per method
% Below zero LN is preferred. Dashed lines mark the benchmarks and n*
figure;
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
for i = 1:numMethods
    if ~hasMethodLog(i)
        continue
    end

    % Calculate BIC_LN - BIC_MLN along nSweep
    sweep_LNminusMLN = 2.*nSweep.*gainMLN(i) - extraParams.*log(nSweep);

    % Plot the difference, the zero line and the benchmarks
    nexttile
    hold on
    plot(nSweep, sweep_LNminusMLN, 'k-', 'LineWidth', 1.5)
    yline(0, 'k:')
    xline(numCores, '--', "cores", 'Color', [0.5 0.5 0.5])
    xline(numCpairs(i), '--', "n_{RP}", 'Color', 'b')
    xline(numWP(i), '--', "n_{WP}", 'Color', [0.5 0.5 0.5])

    % Mark Kish n and n* where they exist
    if isfinite(kishN(i))
        xline(kishN(i), '--', "Kish", 'Color', [0 0.6 0])
    end
    if isfinite(crossoverN(i))
        xline(crossoverN(i), '-', "n*", 'Color', 'r')
    end

    set(gca, 'XScale', 'log')
    ylim(figYLim)
    box on
    title(methodNames(i), 'FontWeight', 'normal')
    xlabel("n")
    ylabel("BIC_{LN} - BIC_{MLN}")
end
