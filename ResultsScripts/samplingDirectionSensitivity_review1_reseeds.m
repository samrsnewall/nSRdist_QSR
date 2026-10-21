%% samplingDirectionSensitivity_review1_reseeds - full 10000-run reseeds of RSR500 for the response to Reviewer 1
%
% Runs calcData again once for each seed in seedList with S.numruns runs,
% making the RSRx restrictions in seedMinAgeDiffs under the three sampling
% orders (no rejection or Bchron), on the cores the manuscript fit used.
% Each seed file is fitted with fitData, and the BICeff of the four
% distributions in each seed is compared with the manuscript fit of all
% S.numruns runs.
%
% Cost: about 2-2.5 h per seed for all nine methods x 69 cores at 10000 runs
% (all three RSRx loops are parfor), plus the scenario step and ~10 min of
% fitting per seed
%
% calcData and fitData are run as scripts in this workspace and pick up the
% calcOverride / fitOverride structs set below; every other setting comes
% from the settings blocks inside those two scripts.
%
% The script cds to the repository root because calcData and fitData use
% paths relative to it.

cd(fileparts(fileparts(mfilename('fullpath'))))
addpath('Functions')

%% Settings
runTag            = "QSR";            % Run the seeds are replicates of (used in the file names)
seedList          = [2 3 4];          % rng seeds of the reseeded calcData runs (the main run used seed 1)
numruns           = 10000;            % S.numruns: Monte Carlo iterations per seed, as in the main run
seedMinAgeDiffs   = [0 500 1000];     % S.rsrMinAgeDiffs: the RSRx restrictions made, under every order in
                                      %   calcData's S.sampleOrders
seedFolder        = "Results/review1_seeds"; % Where the seed calcData files and their fits are saved
pooledFitFolder   = "Results/results_QSR"; % Folder of the manuscript fit of runTag (the pooled
                                      %   fit of all S.numruns runs); its cores are the ones reseeded
weighting         = "depth";  % fitS.weighting: "depth", "age" or "none"; the fit file suffix follows from it
maxFitPoints      = 5e6;      % fitS.maxFitPoints: cap on weighted points passed to the fitters (as in the
                              %   manuscript fit)
weightingPipeline = "sample"; % fitS.weightingPipeline: "replicate" or "sample" (same distribution, much faster)
otherOrderTags    = ["d2s", "random"]; % Suffixes of the orders not used in the manuscript; the manuscript's
                              %   shallow-to-deep methods (no suffix) are marked in the printed summary
tieFactor         = 2;        % Two families count as tied when the gap between their pooled BICeffs is less
                              %   than tieFactor times the sd of that gap across the seeds
fitsOnly          = false;    % If true, calcData is not run: the existing seed files in seedFolder are fitted,
                              %   and any existing seed fit files are deleted first
plotsOnly         = true;    % If true, nothing is made or fitted: the seed fit files already in seedFolder
                              %   are loaded, summarised and plotted
fig2MethodGrid    = ["RSR0",    "RSR0d2s",    "RSR0random"; ...
                     "RSR500",  "RSR500d2s",  "RSR500random"; ...
                     "RSR1000", "RSR1000d2s", "RSR1000random"]; % Methods shown in Fig 2, laid out as this
                              %   grid (rows share a y-axis)
panelLetters      = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]; % Panel labels of Fig 2, row by row

% Set up the file names (the calcData file names without .mat, as fitData
% takes them); fitSuffix must match the rule in fitData.m
fitSuffix     = "_fit_" + weighting + "weight";
pooledFitFile = fullfile(pooledFitFolder, "dataT_" + runTag + fitSuffix + ".mat");
numSeeds      = numel(seedList);
seedFiles     = strings(numSeeds, 1);
seedFitFiles  = strings(numSeeds, 1);
for iSeed = 1:numSeeds
    seedFiles(iSeed)    = fullfile(seedFolder, "dataT_" + runTag + "_seed" + seedList(iSeed));
    seedFitFiles(iSeed) = seedFiles(iSeed) + fitSuffix + ".mat";
end

% calcData saves to Results/<S.dataOutputFile>, so the names given to
% calcData are the seed file names relative to Results/
seedOutputFiles = extractAfter(seedFiles, "Results/");

% calcData and fitData define many short variable names, so the workspace is
% cleared of everything except these between runs
keepVars = {'runTag' 'seedList' 'numruns' 'seedMinAgeDiffs' 'seedFolder' 'weighting' 'maxFitPoints' ...
            'weightingPipeline' 'otherOrderTags' 'tieFactor' 'fitsOnly' 'plotsOnly' 'fig2MethodGrid' ...
            'panelLetters' 'fitSuffix' 'pooledFitFile' 'numSeeds' 'seedFiles' 'seedFitFiles' ...
            'seedOutputFiles' 'fittedCoreNames' 'iSeed' 'keepVars'};

%% Check the seed fit files: they must exist to plot only, and must not exist to fit
for iSeed = 1:numSeeds
    if plotsOnly && ~isfile(seedFitFiles(iSeed))
        error("Fit file not found: " + seedFitFiles(iSeed) + " - set plotsOnly = false to make it")
    elseif ~plotsOnly && isfile(seedFitFiles(iSeed)) && fitsOnly
        delete(seedFitFiles(iSeed))
    elseif ~plotsOnly && isfile(seedFitFiles(iSeed))
        error("Fit file already exists: " + seedFitFiles(iSeed) + " - change seedFolder or delete the file, or set plotsOnly = true to plot it")
    end
end
if ~isfolder(seedFolder)
    mkdir(seedFolder)
end

%% Get the cores the manuscript fit used
% Only these are reseeded, so calcData doesn't sample the other cores in
% the full table, which fitData would leave out anyway (S.coreSelection
% takes a string array of core names)
fittedCoreNames = strings(1, 0);
if ~fitsOnly && ~plotsOnly
    disp("Loading " + pooledFitFile)
    pooledFit       = load(pooledFitFile, "d");
    fittedCoreNames = pooledFit.d.S1.chosenCoreNames;
    clear pooledFit
    disp("Reseeding " + numel(fittedCoreNames) + " cores")
end

%% Run calcData for each seed: RSR500 under every sampling order, nothing else
% oneCoreScenarios still runs for every core (it makes the age-depth
% scenarios RSR500 samples from), and calcData still reaches its final
% addvars and save with S.useBchron false, because the Bchron cells are
% initialised whether or not they are filled
if ~fitsOnly && ~plotsOnly
    for iSeed = 1:numSeeds
        clearvars('-except', keepVars{:})
        seed = seedList(iSeed);

        % Skip a seed whose calcData file already exists (calcData would overwrite it)
        if isfile(seedFiles(iSeed) + ".mat")
            disp("calcData file for seed " + seed + " already exists, not remaking it: " + seedFiles(iSeed))
            continue
        end

        calcOverride.dataOutputFile = seedOutputFiles(iSeed);
        calcOverride.coreSelection  = fittedCoreNames;
        calcOverride.maxNumCores    = [];
        calcOverride.numruns        = numruns;
        calcOverride.rsrMinAgeDiffs = seedMinAgeDiffs;
        calcOverride.useRejection   = false;
        calcOverride.useBchron      = false;

        rng(seed)
        disp("===== calcData (seed " + seed + ") started " + string(datetime("now")))
        calcData
        disp("===== calcData (seed " + seed + ") finished " + string(datetime("now")))
    end
    clearvars('-except', keepVars{:})
end

%% Fit each seed on every core it has, fitting only its RSR columns
% The seed files have no rejection or Bchron columns, so every method list
% is empty and fitData fits just the RSR columns it finds
if ~plotsOnly
    for iSeed = 1:numSeeds
        clearvars('-except', keepVars{:})
        fitOverride.filepath           = seedFiles(iSeed);
        fitOverride.saveFolder         = seedFolder;
        fitOverride.weighting          = weighting;
        fitOverride.maxFitPoints       = maxFitPoints;
        fitOverride.weightingPipeline  = weightingPipeline;
        fitOverride.commonCoresMethods = strings(1, 0);
        fitOverride.lin2014Methods     = strings(1, 0);
        fitOverride.unthinnedMethods   = strings(1, 0);

        disp("===== fitData (seed " + seedList(iSeed) + ") started " + string(datetime("now")))
        fitData
        disp("===== fitData (seed " + seedList(iSeed) + ") finished " + string(datetime("now")))
    end
    clearvars('-except', keepVars{:})
end

%% Get the BICeff of every method and family from the pooled fit and the seed fits
% The pooled fit of all S.numruns runs is the reference; the seeds are the
% replicates. Each fit file is loaded in turn and cleared
distNames = ["LN", "MLN", "Gam", "invGam"];   % Field names in the fit structs
distNamesPlotting = ["LN", "MLN", "G", "IG"]; % dist names for use when plotting...

% Find the methods in the seed fits (the RSR500 columns), leaving out the
% fields that aren't fits
seedFit     = load(seedFitFiles(1), "d");
methodNames = string(fieldnames(seedFit.d.S1));
fitLog      = ~ismember(methodNames, ["fitS", "chooseLog", "numCores", "chosenCoreNames"]) & ~endsWith(methodNames, "IR");
methodNames = methodNames(fitLog);
numMethods  = numel(methodNames);
numDists    = numel(distNames);     
clear seedFit

% Initialise the BICeff arrays
replicateBICeff = NaN(numMethods, numDists, numSeeds);
referenceBICeff = NaN(numMethods, numDists);

% Get the BICeff of each seed
for iSeed = 1:numSeeds
    disp("Loading " + seedFitFiles(iSeed))
    seedFit = load(seedFitFiles(iSeed), "d");
    for iMethod = 1:numMethods
        for iDist = 1:numDists
            methodFit = seedFit.d.S1.(methodNames(iMethod));
            replicateBICeff(iMethod, iDist, iSeed) = methodFit.(distNames(iDist)).fitInfo.BICtaeheefix;
        end
    end
    clear seedFit
end

% Get the BICeff of the pooled fit
disp("Loading " + pooledFitFile)
pooledFit = load(pooledFitFile, "d");
for iMethod = 1:numMethods
    if ~isfield(pooledFit.d.S1, methodNames(iMethod))
        warning(methodNames(iMethod) + " is not in the pooled fit file")
        continue
    end
    for iDist = 1:numDists
        methodFit = pooledFit.d.S1.(methodNames(iMethod));
        referenceBICeff(iMethod, iDist) = methodFit.(distNames(iDist)).fitInfo.BICtaeheefix;
    end
end
clear pooledFit

%% Summarise the BICeff of each method across the replicates
% For each family: the reference BICeff and the median, min, max and sd
% across the replicates. Then the reference ranking (lowest BICeff first)
% and the fraction of replicates that rank all four families the same way.
% Then, for the LN-MLN and Gam-invGam pairs, the gap in the reference fit,
% the mean and sd of the gap across the replicates, and whether the pair is
% tied (|reference gap| < tieFactor x replicate sd). With only numSeeds
% replicates the sd, and so the tied flag, is rough
numReplicates = size(replicateBICeff, 3);
gapPairs      = ["LN", "MLN"; "Gam", "invGam"];   % Each row is one pair; the gap is first minus second
BICeffT       = table(methodNames, 'VariableNames', "method");

% Get the reference value and the spread across the replicates of each family
for iDist = 1:numDists
    distBICeff = squeeze(replicateBICeff(:, iDist, :));
    BICeffT.(distNames(iDist) + "_reference")       = referenceBICeff(:, iDist);
    BICeffT.(distNames(iDist) + "_replicateMedian") = median(distBICeff, 2);
    BICeffT.(distNames(iDist) + "_replicateMin")    = min(distBICeff, [], 2);
    BICeffT.(distNames(iDist) + "_replicateMax")    = max(distBICeff, [], 2);
    BICeffT.(distNames(iDist) + "_replicateSD")     = std(distBICeff, 0, 2);
end

% Get the reference ranking of each method and how many replicates agree with it
referenceRanking    = strings(numMethods, 1);
fracReplicatesAgree = NaN(numMethods, 1);
for iMethod = 1:numMethods
    [~, referenceOrder]       = sort(referenceBICeff(iMethod, :));
    referenceRanking(iMethod) = strjoin(distNames(referenceOrder), " > ");

    % Count the replicates whose ranking is the same as the reference ranking
    numAgree = 0;
    for iRep = 1:numReplicates
        [~, replicateOrder] = sort(replicateBICeff(iMethod, :, iRep));
        if isequal(replicateOrder, referenceOrder)
            numAgree = numAgree + 1;
        end
    end
    fracReplicatesAgree(iMethod) = numAgree./numReplicates;
end
BICeffT.referenceRanking    = referenceRanking;
BICeffT.fracReplicatesAgree = fracReplicatesAgree;

% Get the gap of each pair, and whether the pair is tied
tiedText = strings(numMethods, 1);
for iPair = 1:size(gapPairs, 1)
    % Find the two families of this pair
    firstIdx  = find(distNames == gapPairs(iPair, 1));
    secondIdx = find(distNames == gapPairs(iPair, 2));
    pairName  = gapPairs(iPair, 1) + "_" + gapPairs(iPair, 2);

    % Calculate the gaps, in the reference fit and in each replicate
    referenceGap = referenceBICeff(:, firstIdx) - referenceBICeff(:, secondIdx);
    replicateGap = squeeze(replicateBICeff(:, firstIdx, :) - replicateBICeff(:, secondIdx, :));
    replicateGapSD = std(replicateGap, 0, 2);
    tiedLog        = abs(referenceGap) < tieFactor.*replicateGapSD;

    % Store them, and note the tied pairs for the printed summary
    BICeffT.(pairName + "_referenceGap")  = referenceGap;
    BICeffT.(pairName + "_replicateMean") = mean(replicateGap, 2);
    BICeffT.(pairName + "_replicateSD")   = replicateGapSD;
    BICeffT.(pairName + "_tied")          = tiedLog;
    tiedText(tiedLog) = tiedText(tiedLog) + " " + gapPairs(iPair, 1) + "-" + gapPairs(iPair, 2);
end

% Save the summary
writetable(BICeffT, fullfile(seedFolder, "BICeff_seeds_" + runTag + ".csv"))

% Print one line per method: the reference ranking, the fraction of replicates
% that agree with it, and which pairs are tied. Methods in the manuscript's
% order (shallow to deep, no order suffix) are marked with a *
disp("BICeff ranking of the pooled fit (best first), fraction of the " + numReplicates + ...
    " seeds with the same ranking, and the tied pairs (tieFactor = " + tieFactor + "):")
for iMethod = 1:numMethods
    if tiedText(iMethod) == ""
        tiedText(iMethod) = " none";
    end
    chosenMark = " ";
    if ~endsWith(methodNames(iMethod), otherOrderTags)
        chosenMark = "*";
    end
    disp(" " + chosenMark + " " + methodNames(iMethod) + ": " + referenceRanking(iMethod) + "   agree " + ...
        compose("%.2f", fracReplicatesAgree(iMethod)) + "   tied:" + tiedText(iMethod))
end

%% Fig 2: BICeff of each family across the seeds, one panel per method laid out as fig2MethodGrid
% Rows are the RSRx restrictions and columns the sampling orders; the panels
% of a row share a y-axis (the restrictions sit at very different BICeff
% levels). The manuscript run (seed 1, the pooled fit) is one more
% replicate of the same size, so it is drawn as a dot like the seeds (an open
% red circle) rather than as a reference, and the short line is the mean of
% all the replicates including it
[numGridRows, numGridCols] = size(fig2MethodGrid);
allSeedBICeff  = cat(3, replicateBICeff, referenceBICeff);   % numMethods x numDists x (numSeeds + 1)
meanSeedBICeff = mean(allSeedBICeff, 3);
figure('Units', 'centimeters', 'Position', [2 2 6.*numGridCols 6.*numGridRows]);
tiledlayout(numGridRows, numGridCols, 'TileSpacing', 'compact', 'Padding', 'compact')
fig2Axes = gobjects(numGridRows, numGridCols);
for iRow = 1:numGridRows
    for iCol = 1:numGridCols
        % Find this method in the seed fits
        gridMethodName = fig2MethodGrid(iRow, iCol);
        iMethod        = find(methodNames == gridMethodName);
        iPanel         = (iRow-1).*numGridCols + iCol;
        if isempty(iMethod)
            warning(gridMethodName + " is not in the seed fits, so its Fig 2 panel is empty")
            nexttile(iPanel)
            continue
        end

        fig2Axes(iRow, iCol) = nexttile(iPanel);
        hold on
        for iDist = 1:numDists
            % Plot the reseeds, the manuscript run and the mean of all of them
            distBICeff = squeeze(replicateBICeff(iMethod, iDist, :));
            plot(iDist, referenceBICeff(iMethod, iDist), 'o', 'Color', 'r', 'MarkerSize', 4, 'LineWidth', 1)
            plot(iDist.*ones(numReplicates, 1), distBICeff, 'x','Color', 'k', 'MarkerSize', 8)
            plot(iDist + [-0.3 0.3], meanSeedBICeff(iMethod, iDist).*[1 1], 'k-', 'LineWidth', 1.5)
        end
        % Label the panel inside the axes ("A: RSR0") instead of a title
        text(0.04, 0.92, panelLetters(iPanel) + ": " + methodNames(iMethod), 'Units', 'normalized', ...
            'FontSize', 9, 'FontWeight', 'bold')
        box on
        xlim([0.5 numDists + 0.5])
        xticks(1:numDists)
        xticklabels(distNamesPlotting)
        if iCol == 1
            ylabel("BICeff")
        end
        if iRow < numGridRows
            set(gca, 'XTickLabel', [])
        end
    end

    % Share the y-axis along this row
    rowAxes = fig2Axes(iRow, :);
    linkaxes(rowAxes(isgraphics(rowAxes)), 'y')
end

% Put one horizontal legend along the bottom of the figure, below the x tick
% labels of the bottom row, from the series of the first panel
firstAxes  = fig2Axes(isgraphics(fig2Axes));
seedLegend = legend(firstAxes(1), ["manuscript run (seed 1)", "reseeds", "mean"], ...
    'Orientation', 'horizontal', 'Box', 'on', 'FontSize', 8);
seedLegend.Layout.Tile = 'south';
