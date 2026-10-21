%% results_QSR - run the nSRdist pipeline and make the manuscript figures
% This script runs the whole nSRdist pipeline and then makes the manuscript
% figures and tables. It runs calcData on the production core set, then
% fitData once for each weighting in runS.weightings (the first one also
% gets the individual-run fits), then plots. With runS.plotsOnly = true,
% calcData and fitData are skipped and the saved files are loaded and
% plotted.
%
% calcData and fitData are run as scripts in this workspace and pick up the
% calcOverride / fitOverride structs set below; every other setting comes
% from the settings blocks inside calcData and fitData. Their workspaces are
% cleared between runs, keeping only runS and keepVars.
%
% The manuscript methods are stored in the fit file as:
%   BMedian -> BMedianLin, BSamp -> BSampLin (both with the Lin2014 age
%   filter), RSR0/RSR500/RSR1000 (shallow-to-deep sampling), and their
%   individual-run fits BSampLinIR, RSR500IR, ...
%
% Histograms of pooled fits are scaled to cm of sediment per run (counts x
% sedLength / number of points), so the y axes stay in cm whether or not
% fitData capped the number of points (runS.maxFitPoints). Individual-run
% histograms (lnSRHistCounts) are already in cm.
%
% The script cds to the repository root because calcData and fitData use
% paths relative to it. Ages in kyr unless stated; depths in cm.

cd(fileparts(fileparts(mfilename('fullpath'))))
addpath("Functions")

%% Settings
runS.runTag            = "QSR";           % Used in the calcData and fit file names
runS.coreSelection     = "good";          % S.coreSelection: "good" for the full production core set, or a string array of core names
runS.maxNumCores       = [];              % S.maxNumCores: keep only the first this-many chosen cores alphabetically ([] = all)
runS.numruns           = 10000;           % S.numruns: Monte Carlo iterations for BSamp, RSRx and rejection

runS.calcSeed          = 1;               % rng seed before calcData
runS.weightings        = ["depth", "age", "none"];   % fitS.weighting for each fitData run; the first is the main fit
                                          %   (all methods and individual-run fits), the others fit only the
                                          %   RSR columns for the weighting comparison
runS.maxFitPoints      = 5e6;             % fitS.maxFitPoints: cap on weighted points passed to the fitters
                                          %   (with all cores at 10000 runs an uncapped pooled fit takes hours)
runS.weightingPipeline = "sample";        % fitS.weightingPipeline: "replicate" or "sample" (same distribution, much faster)
runS.numIRruns         = 1000;            % fitS.OneRun.numruns: individual runs fitted in the main fit
runS.fitFolder         = "Results/results_QSR";   % Where this script's fit files are saved (fitOverride.saveFolder)
runS.plotsOnly         = true;           % If true, calcData and fitData are skipped and the saved files are plotted
runS.reuseCalcFile     = true;            % If true and the calcData file already exists (e.g. a run that stopped in
                                          %   fitData), calcData is skipped and fitData uses that file
% Set up the file names; the fit suffix must match the rule in fitData.m
runS.calcFile     = "Results/dataT_" + runS.runTag;
runS.coreInfoFile = "Results/dataT_" + runS.runTag + "_coreInfo.mat";
runS.fitFiles     = fullfile(runS.fitFolder, "dataT_" + runS.runTag + "_fit_" + runS.weightings + "weight.mat");

% calcData and fitData define many short variable names, so the workspace is
% cleared of everything except these between runs
keepVars = {'runS', 'keepVars'};

%% Check the files: they must exist to plot only, and must not exist to run
if runS.plotsOnly
    if ~isfile(runS.calcFile + ".mat") || any(~isfile(runS.fitFiles))
        error("Missing calcData or fit file for " + runS.runTag + " - set runS.plotsOnly = false to make them")
    end
else
    if isfile(runS.calcFile + ".mat") && ~runS.reuseCalcFile
        error("calcData file already exists: " + runS.calcFile + " - change runS.runTag, or set runS.reuseCalcFile = true")
    end
    if any(isfile(runS.fitFiles))
        error("A fit file already exists in " + runS.fitFolder + " - change runS.runTag or runS.fitFolder")
    end
end
if ~isfolder(runS.fitFolder)
    mkdir(runS.fitFolder)
end

%% Run calcData on the production core set
% Skipped if runS.reuseCalcFile and the calcData file already exists
runS.runCalcData = ~runS.plotsOnly && ~(runS.reuseCalcFile && isfile(runS.calcFile + ".mat"));
if ~runS.plotsOnly && ~runS.runCalcData
    disp("Reusing existing calcData file " + runS.calcFile)
end
if runS.runCalcData
    clearvars('-except', keepVars{:})
    calcOverride.dataOutputFile = "dataT_" + runS.runTag;
    calcOverride.coreSelection  = runS.coreSelection;
    calcOverride.maxNumCores    = runS.maxNumCores;
    calcOverride.numruns        = runS.numruns;

    rng(runS.calcSeed)
    disp("===== calcData started " + string(datetime("now")))
    calcData
    disp("===== calcData finished " + string(datetime("now")))

    % Save the core metadata the figures need in a small file, so plotting
    % doesn't have to load dataT (over 1 GB)
    coreInfoT = dataT(:, ["cores", "lats", "longs", "depths", "meanSR"]);
    save(runS.coreInfoFile, "coreInfoT")
    clearvars('-except', keepVars{:})
end

%% Run fitData for each weighting
if ~runS.plotsOnly
    for iWeight = 1:numel(runS.weightings)
        loopKeepVars = [keepVars, {'iWeight'}];
        clearvars('-except', loopKeepVars{:})
        fitOverride.filepath          = runS.calcFile;
        fitOverride.saveFolder        = runS.fitFolder;
        fitOverride.weighting         = runS.weightings(iWeight);
        fitOverride.maxFitPoints      = runS.maxFitPoints;
        fitOverride.weightingPipeline = runS.weightingPipeline;

        % The main fit (first weighting) fits every method with individual
        % runs; the others fit only the RSR columns (fitData always fits
        % every RSR column), which is all the weighting comparison uses
        if iWeight == 1
            fitOverride.fitIRs        = true;
            fitOverride.OneRunNumruns = runS.numIRruns;
        else
            fitOverride.fitIRs           = false;
            fitOverride.lin2014Methods   = strings(1, 0);
            fitOverride.unthinnedMethods = strings(1, 0);
        end

        disp("===== fitData (" + runS.weightings(iWeight) + " weighting) started " + string(datetime("now")))
        fitData
        disp("===== fitData (" + runS.weightings(iWeight) + " weighting) finished " + string(datetime("now")))
    end
    clearvars('-except', keepVars{:})
end

%% Load the fit results and the core metadata
% Load the main fit (first weighting)
mainFit = load(runS.fitFiles(1), "d");
d       = mainFit.d;

% Load the core metadata, making the small file from dataT if it is missing
% (e.g. if calcData was run outside this script)
if ~isfile(runS.coreInfoFile)
    disp("Making " + runS.coreInfoFile + " from " + runS.calcFile + " (loads dataT, slow)")
    calcResults = load(runS.calcFile, "dataT");
    coreInfoT   = calcResults.dataT(:, ["cores", "lats", "longs", "depths", "meanSR"]);
    save(runS.coreInfoFile, "coreInfoT")
    clear calcResults
end
load(runS.coreInfoFile, "coreInfoT")

% Keep the metadata of the fitted cores only
fittedLog = ismember(string(coreInfoT.cores), d.S1.chosenCoreNames);
dataTbl   = coreInfoT(fittedLog, :);
disp(num2str(height(dataTbl)) + " cores in the main fit")

% Set up the names of the manuscript methods in this fit file
bmedianName = "BMedianLin";
bsampName   = "BSampLin";
rsr0Name    = "RSR0";
rsr500Name  = "RSR500";
rsr1000Name = "RSR1000";

%% Set text size for figures
set(groot, 'DefaultAxesFontSize', 12)
set(groot, 'DefaultTextFontSize', 12)

%% Set up some general histogram settings
bw            = 0.1;           % Bin width in log(nSR) space
logBinEdges   = -5:bw:5;
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end)).*0.5;
dtBinEdges    = 0:0.1:10;      % kyr; matches the 100 yr bins of agediffsWC / agediffsWbC

% Set up the colours and line styles of the four fitted distributions,
% from the colour-blind friendly six-colour scheme of Petroff (2021,
% arXiv:2107.02270). The lognormal family is warm (LN red, MLN orange) and
% the gamma family cool (G purple, IG blue), with the colours ordered by
% how visible they are on the grey histograms: LN matters most, G least,
% and G and IG rarely overlap so they can be close
fitColours.LN     = [228  37  54]./255;   % Red
fitColours.MLN    = [248 156  32]./255;   % Orange
fitColours.Gam    = [150  74 139]./255;   % Purple
fitColours.invGam = [ 87 144 252]./255;   % Blue
fitLineStyles.LN     = '-';
fitLineStyles.MLN    = '-';
fitLineStyles.Gam    = '-';
fitLineStyles.invGam = '-';

%% Input and analyse BIGMACS data
% Load BIGMACS files
BMlognorm    = readtable("BIGMACSdata/lognormal.txt");
BM.nSR.x     = BMlognorm.Var1';
BM.nSR.px    = BMlognorm.Var2;
BM.weightedC = readmatrix("BIGMACSdata/Lin2014_sedrateratio_cm_wo_NaN.txt"); % all weighted counts
BM.weightedC = BM.weightedC(:,4);
BM.sedLength   = length(BM.weightedC); % cm; every nSR value represents 1 cm of sediment
BM.sedTimeSpan = 544;                  % kyr; from Lin et al. 2014
BM.numCpairs   = 343;                  % known from running code used for Lin2014 (not available to others)

% Convert the BIGMACS nSR pdf to log(nSR)
[BM.lnSR.x, BM.lnSR.px] = px_to_pfx(BM.nSR.x, BM.nSR.px, @log);

%% Figure 1: map of cores, histogram of water depths, histogram of mean SRs
% Define tick locations
lonTicks = -180:20:180;
latTicks = -40:20:40;

figure;

% Map panel
subplot(2,2,[1 2])
ax = axesm('mercator', 'MapLatLimit', [-45 45], 'MapLonLimit', [-180 180], ...
    'MeridianLabel', 'on', 'ParallelLabel', 'on', 'Frame', 'off', ...
    'LabelFormat', 'none', 'MLabelLocation', 20, 'PLabelLocation', 20, ...
    'MLabelParallel', -45);
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8])
load coastlines coastlat coastlon
plotm(coastlat, coastlon, 'Color', 'k')
tightmap
plotm(dataTbl.lats, dataTbl.longs, 'ks', 'LineWidth', 1)
set(gca, 'OuterPosition', [0 0.4 1 0.6])  % [left, bottom, width, height]
text(0.02, 0.9, 'A', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% Draw on ticks
xl       = xlim;
yl       = ylim;
xTickLen = 0.05;
yTickLen = 0.05;
proj     = gcm;
[lonTick_x, ~] = projfwd(proj, repmat(latTicks(1), size(lonTicks)), lonTicks);
[~, latTick_y] = projfwd(proj, latTicks, repmat(lonTicks(1), size(latTicks)));

% Draw longitude ticks at bottom and top
for i = 1:length(lonTicks)
    line([lonTick_x(i) lonTick_x(i)], [yl(1) yl(1)+yTickLen], 'Color', 'k', 'LineWidth', 0.5)
    line([lonTick_x(i) lonTick_x(i)], [yl(2)-yTickLen yl(2)], 'Color', 'k', 'LineWidth', 0.5)
end

% Draw latitude ticks at left and right
for i = 1:length(latTicks)
    line([xl(1) xl(1)+xTickLen], [latTick_y(i) latTick_y(i)], 'Color', 'k', 'LineWidth', 0.5)
    line([xl(2)-xTickLen xl(2)], [latTick_y(i) latTick_y(i)], 'Color', 'k', 'LineWidth', 0.5)
end

% Rotate longitude labels and shift them down
xticklabelHandles = findall(gcf, 'Type', 'text', 'Tag', 'MLabel');
numXLabels        = numel(xticklabelHandles);
xtickPosition     = cell2mat(get(xticklabelHandles, 'Position'));
hold(ax, 'on')
xTickLabels = get(xticklabelHandles, 'String');
xTickLabels = cellfun(@(c)c(2), xTickLabels);
set(xticklabelHandles, {'String'}, xTickLabels)
set(xticklabelHandles, 'Rotation', 90)
set(xticklabelHandles, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'Middle')
xtickPosition(:,2) = xtickPosition(:,2) - range(ax.YLim)*.03;
set(xticklabelHandles, {'Position'}, mat2cell(xtickPosition, ones(numXLabels,1), 3))
tightmap

% Shift latitude labels right
yticklabelHandles = findall(gcf, 'Type', 'text', 'Tag', 'PLabel');
numYLabels        = numel(yticklabelHandles);
ytickPosition     = cell2mat(get(yticklabelHandles, 'Position'));
ytickPosition(:,1) = ytickPosition(:,1) + range(ax.XLim)*.02;
set(yticklabelHandles, {'Position'}, mat2cell(ytickPosition, ones(numYLabels,1), 3))
tightmap

% Histogram of core water depths
subplot(2,2,3)
histogram(dataTbl.depths./1000, 0:0.5:6, 'FaceColor', [0.8 0.8 0.8])
xlim([0 6])
xlabel("Depth (km)")
ylabel("Counts")
xticks(0:1:6)
ylims = ylim;
ylim(ylims.*1.1)
set(gca, 'Position', [0.13 0.25 0.33 0.25])
text(0.02, 0.9, 'B', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% Histogram of core average SRs
subplot(2,2,4)
histogram(dataTbl.meanSR, 0:5:90, 'FaceColor', [0.8 0.8 0.8])
xlabel('Average SR (cm/kyr)')
ylabel('Counts')
xlim([0 90])
ylims = ylim;
ylim(ylims.*1.1)
set(gca, 'Position', [0.57 0.25 0.33 0.25])
text(0.02, 0.9, 'C', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

%% Create table of useful information for each method
dStrus        = {BM, d.S1.(bmedianName), d.S1.(bsampName), d.S1.(rsr0Name), d.S1.(rsr500Name), d.S1.(rsr1000Name)};
dStrusStrings = ["BIGMACS", bmedianName, bsampName, rsr0Name, rsr500Name, rsr1000Name];

% Initialise the table columns
MeanAgePairsT = NaN(length(dStrus),1);
MeanSedLength = NaN(length(dStrus),1);
MeanSedTime   = NaN(length(dStrus),1);
nSR_median    = NaN(length(dStrus),1);
nSR_min       = NaN(length(dStrus),1);
nSR_max       = NaN(length(dStrus),1);
nSR_95lo      = NaN(length(dStrus),1);
nSR_95hi      = NaN(length(dStrus),1);

% Fill them in for each method. BIGMACS sedTimeSpan is in kyr; the fits'
% sedTimeSpan is in yr (sum of dt_yr per run), so convert those to kyr
for i = 1:length(dStrus)
    dStru = dStrus{i};
    allWC = dStru.weightedC;
    MeanSedLength(i,1) = dStru.sedLength;
    if i == 1
        MeanSedTime(i,1) = dStru.sedTimeSpan;
    else
        MeanSedTime(i,1) = dStru.sedTimeSpan./1000;
    end
    MeanAgePairsT(i,1) = mean(dStru.numCpairs);

    % Get the nSR summary statistics
    allWCsorted     = sort(allWC);
    numWC           = length(allWC);
    nSR_median(i,1) = median(allWC);
    nSR_min(i,1)    = min(allWC);
    nSR_max(i,1)    = max(allWC);
    nSR_95lo(i,1)   = allWCsorted(ceil(0.025*numWC));
    nSR_95hi(i,1)   = allWCsorted(ceil(0.975*numWC));
end
dA.summT = table(MeanAgePairsT, MeanSedLength, MeanSedTime, nSR_median, ...
    nSR_min, nSR_max, nSR_95lo, nSR_95hi, 'RowNames', dStrusStrings);

%% Table 1: amount of data/sediment used in each method (sediment time span in kyr)
dA.summT(:, 1:3)

%% Figure 2: BIGMACS and BMedian nSR with the BIGMACS pdf, and histograms of dt
% Load the dts of the BIGMACS data and get their weighted histogram
load("Results/Lin2014_dts.mat", "asr")
BM_dt_weightedBC = makeWeightedBinCounts(asr(:,2), asr(:,3), dtBinEdges);

% Set up the histogram colours and bar edge widths: the retained data
% (BIGMACS, BMedianLin) in dark grey with thicker edges, the unfiltered
% BMedian behind BMedianLin in light grey with thin edges
histColour          = [0.75 0.75 0.75]; % Face colour of the retained histograms
unfilteredColour    = [0.95 0.95 0.95]; % Face colour of the unfiltered histograms behind BMedianLin
histLineWidth       = 0.5;             % Bar edge width of the retained histograms
unfilteredLineWidth = 0.5;           % Bar edge width of the unfiltered histograms behind them
unfilteredEdgeColour = [0.5 0.5 0.5]; % Bar edge colour of the unfiltered histograms
unfilteredLineStyle  = '-';          % Bar edge line style of the unfiltered histograms

% Set up figure
fA = figure;
fA.Position = [100, 100, 1000, 400]; % [left, bottom, width, height]

% Plot histogram of BIGMACS weighted nSR counts (each value is 1 cm)
subplot(2,2,1)
BM.logHC = histcounts(log(BM.weightedC), logBinEdges);
normF    = 1./sum(BM.logHC.*bw);
hold on
histogram('BinCounts', BM.logHC, 'BinEdges', logBinEdges, 'FaceColor', histColour, 'FaceAlpha', 1, 'LineWidth', histLineWidth)
plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
text(0.02, 0.9, 'A: BIGMACS', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
xlim([-4 4])
ylabel("cm")
xlabel("log(NSR)")

% Plot histogram of BIGMACS weighted dts
subplot(2,2,2)
hold on
histogram('BinCounts', BM_dt_weightedBC, 'BinEdges', dtBinEdges, 'FaceColor', histColour, 'FaceAlpha', 1, 'LineWidth', histLineWidth)
text(0.02, 0.9, 'B: BIGMACS', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
ylabel("cm")
xlabel("\Deltat (kyr)")
xlim([0 6])

% Get the unfiltered BMedian histograms (no Lin2014 dt filter; d.S1.BMedian
% is made by having BMedian in fitS.unthinnedMethods), scaled to cm per run
% the same way, to plot behind the BMedianLin ones so the trimmed data can
% be seen
bmedianUnfiltered  = d.S1.BMedian;
unfiltered_logHC   = histcounts(log(bmedianUnfiltered.weightedC), logBinEdges);
unfiltered_logHC   = unfiltered_logHC.*bmedianUnfiltered.sedLength./numel(bmedianUnfiltered.weightedC);
unfiltered_dtHC    = bmedianUnfiltered.agediffsWC.*bmedianUnfiltered.sedLength./sum(bmedianUnfiltered.agediffsWC);

% Plot histogram of BMedian weighted nSR counts, scaled to cm per run, with
% the unfiltered BMedian behind it
subplot(2,2,3)
bmedianFit    = d.S1.(bmedianName);
BMedian_logHC = histcounts(log(bmedianFit.weightedC), logBinEdges);
BMedian_logHC = BMedian_logHC.*bmedianFit.sedLength./numel(bmedianFit.weightedC);
normF         = 1./sum(BMedian_logHC.*bw);
hold on
histogram('BinCounts', unfiltered_logHC, 'BinEdges', logBinEdges, 'FaceColor', unfilteredColour, 'FaceAlpha', 1, 'LineWidth', unfilteredLineWidth, ...
        'EdgeColor', unfilteredEdgeColour, 'LineStyle', unfilteredLineStyle)
histogram('BinCounts', BMedian_logHC, 'BinEdges', logBinEdges, 'FaceColor', histColour, 'FaceAlpha', 1, 'LineWidth', histLineWidth);
plot(BM.lnSR.x, BM.lnSR.px*(1/normF), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')
text(0.02, 0.9, 'C: BMedian', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
xlim([-4 4])
ylabel("cm")
xlabel("log(NSR)")

% Plot histogram of BMedian weighted dts, with the unfiltered BMedian
% behind it
subplot(2,2,4)
hold on
histogram('BinCounts', unfiltered_dtHC, 'BinEdges', dtBinEdges, 'FaceColor', unfilteredColour, 'FaceAlpha', 1, 'LineWidth', unfilteredLineWidth, ...
        'EdgeColor', unfilteredEdgeColour, 'LineStyle', unfilteredLineStyle)
histogram('BinCounts', bmedianFit.agediffsWC, 'BinEdges', dtBinEdges, 'FaceColor', histColour, 'FaceAlpha', 1, 'LineWidth', histLineWidth)
text(0.02, 0.9, 'D: BMedian', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
box on
ylabel("cm")
xlabel("\Deltat (kyr)")
xlim([0 6])

% Set the same y-axis limits on the two nSR panels and on the two dt panels
for iCol = 1:2
    subplot(2,2,iCol)
    topYlim = max(ylim);
    subplot(2,2,iCol+2)
    topYlim = max(topYlim, max(ylim));
    subplot(2,2,iCol)
    ylim([0 topYlim])
    subplot(2,2,iCol+2)
    ylim([0 topYlim])
end

%% Figure 3: all sampling approaches, pooled histograms of log(nSR) and dt
% One row per method: BSamp, RSR0, RSR500, RSR1000 (pooled fits). Both
% histograms are scaled to cm per run: the nSR counts by sedLength over the
% number of points the fitters saw (as in Figure 4), and the dt counts by
% sedLength over their total (agediffsWC is summed over all the runs).
% The BSamp row also shows the unfiltered BSamp fit (no Lin2014 dt filter)
% behind the BSampLin histograms, so the trimmed data can be seen, and the
% RSR500 and RSR1000 rows show the outline of the RSR0 histograms, so the
% effect of the minimum dt can be seen
pooledNames    = [bsampName, rsr0Name, rsr500Name, rsr1000Name];
rowNames       = ["BSamp", "RSR0", "RSR500", "RSR1000"]; % Manuscript names for the panel labels
panelLetters   = ["A", "B", "C", "D", "E", "F", "G", "H"];
histColour          = [0.75 0.75 0.75]; % Face colour of the pooled histograms
unfilteredColour    = [0.95 0.95 0.95]; % Face colour of the unfiltered BSamp histogram behind BSampLin
histLineWidth       = 0.5;             % Bar edge width of the pooled histograms
unfilteredLineWidth = 0.5;           % Bar edge width of the unfiltered BSamp histogram
unfilteredEdgeColour = [0.5 0.5 0.5]; % Bar edge colour of the unfiltered BSamp histogram
unfilteredLineStyle  = '-';          % Bar edge line style of the unfiltered BSamp histogram
rsr0OutlineColour    = [0 0 0];      % Line colour of the RSR0 outline on the RSR500 and RSR1000 rows
rsr0OutlineWidth     = 1;            % Line width of the RSR0 outline
bsampUnfiltered  = d.S1.BSamp;
rsr0Fit          = d.S1.(rsr0Name);

% Get the unfiltered BSamp histograms, scaled to cm per run the same way
unfiltered_logHC = histcounts(log(bsampUnfiltered.weightedC), logBinEdges);
unfiltered_logHC = unfiltered_logHC.*bsampUnfiltered.sedLength./numel(bsampUnfiltered.weightedC);
unfiltered_dtHC  = bsampUnfiltered.agediffsWC.*bsampUnfiltered.sedLength./sum(bsampUnfiltered.agediffsWC);

% Get the RSR0 histograms, scaled to cm per run the same way, for the
% outlines on the RSR500 and RSR1000 rows
rsr0_logHC = histcounts(log(rsr0Fit.weightedC), logBinEdges);
rsr0_logHC = rsr0_logHC.*rsr0Fit.sedLength./numel(rsr0Fit.weightedC);
rsr0_dtHC  = rsr0Fit.agediffsWC.*rsr0Fit.sedLength./sum(rsr0Fit.agediffsWC);

fB = figure;
fB.Position = [100, 100, 1000, 800]; % [left, bottom, width, height]
for iMethod = 1:numel(pooledNames)
    pooledFit = d.S1.(pooledNames(iMethod));

    % Plot the pooled log(nSR) histogram, with the unfiltered BSamp behind
    % it on the BSamp row
    subplot(4,2,2*iMethod-1)
    pooled_logHC = histcounts(log(pooledFit.weightedC), logBinEdges);
    pooled_logHC = pooled_logHC.*pooledFit.sedLength./numel(pooledFit.weightedC);
    hold on
    box on
    if iMethod == 1
        histogram('BinCounts', unfiltered_logHC, 'BinEdges', logBinEdges, 'FaceColor', unfilteredColour, 'FaceAlpha', 1, 'LineWidth', unfilteredLineWidth, ...
        'EdgeColor', unfilteredEdgeColour, 'LineStyle', unfilteredLineStyle)
    end
    histogram('BinCounts', pooled_logHC, 'BinEdges', logBinEdges, 'FaceColor', histColour, 'FaceAlpha', 1, 'LineWidth', histLineWidth)
    if iMethod >= 3
        histogram('BinCounts', rsr0_logHC, 'BinEdges', logBinEdges, 'DisplayStyle', 'stairs', ...
            'EdgeColor', rsr0OutlineColour, 'LineWidth', rsr0OutlineWidth)
    end
    xlim([-2.5 2.5])
    xlabel("log(NSR)")
    ylabel("cm per run")
    text(0.02, 0.9, panelLetters(2*iMethod-1) + ": " + rowNames(iMethod), 'Units', 'normalized', ...
        'FontSize', 12, 'FontWeight', 'bold')

    % Plot the pooled dt histogram, with the unfiltered BSamp behind it on
    % the BSamp row
    subplot(4,2,2*iMethod)
    pooled_dtHC = pooledFit.agediffsWC.*pooledFit.sedLength./sum(pooledFit.agediffsWC);
    hold on
    box on
    if iMethod == 1
        histogram('BinCounts', unfiltered_dtHC, 'BinEdges', dtBinEdges, 'FaceColor', unfilteredColour, 'FaceAlpha', 1, 'LineWidth', unfilteredLineWidth, ...
        'EdgeColor', unfilteredEdgeColour, 'LineStyle', unfilteredLineStyle)
    end
    histogram('BinCounts', pooled_dtHC, 'BinEdges', dtBinEdges, 'FaceColor', histColour, 'FaceAlpha', 1, 'LineWidth', histLineWidth)
    if iMethod >= 3
        histogram('BinCounts', rsr0_dtHC, 'BinEdges', dtBinEdges, 'DisplayStyle', 'stairs', ...
            'EdgeColor', rsr0OutlineColour, 'LineWidth', rsr0OutlineWidth)
    end
    xlim([0 6])
    xlabel("\Deltat (kyr)")
    ylabel("cm per run")
    text(0.02, 0.9, panelLetters(2*iMethod) + ": " + rowNames(iMethod), 'Units', 'normalized', ...
        'FontSize', 12, 'FontWeight', 'bold')
end

% Make all the y limits the same within each column of plots, set by the
% tallest panel
fig2Ylims = [4000, 2000];
for iCol = 1:2
    maxYlim = 0;
    for i = iCol:2:8
        subplot(4,2,i)
        %maxYlim = max(maxYlim, max(ylim));
        maxYlim = fig2Ylims(iCol)./1.05;
    end
    for i = iCol:2:8
        subplot(4,2,i)
        ylim([0 maxYlim.*1.05])
    end
end

%% Table 2: BIC of the fitted distributions for each sampling method
bicMethods = [bmedianName, bsampName, rsr0Name, rsr500Name, rsr1000Name];
distNames  = ["LN", "MLN", "Gam", "invGam"];
bicValues  = NaN(numel(bicMethods), numel(distNames));
for i = 1:numel(bicMethods)
    for j = 1:numel(distNames)
        bicValues(i,j) = d.S1.(bicMethods(i)).(distNames(j)).fitInfo.BICtaeheefix;
    end
end
BIC_table = array2table(round(bicValues), 'VariableNames', distNames, 'RowNames', bicMethods)

%% Figure 4: pooled histogram of each method with its four fitted distributions
% One panel per method: BIGMACS, BMedian, BSamp, RSR0, RSR500, RSR1000.
% The BIGMACS data (Lin2014 weighted counts, one per cm) are not fitted in
% fitData, so the four distributions are fitted here with the same fitters
% ARfitdists uses, on the same grid d.x, with BM.numCpairs as the effective
% number of observations. The published BIGMACS lognormal is dashed in
% every panel for reference
fitMethodNames = ["BIGMACS", bmedianName, bsampName, rsr0Name, rsr500Name, rsr1000Name];
panelNames     = ["BIGMACS", "BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"]; % Manuscript names for the panel labels
panelLetters   = ["A", "B", "C", "D", "E", "F"];
fig4Layout     = [2, 3];                   % Panels as [rows, columns]
fig4Position   = [100, 100, 1300, 600];    % [left, bottom, width, height]
fig4Xlim       = [-2.5 2.5];                   % log(nSR) axis limits
logFH          = @(x) log(x);

% Fit the four distributions to the BIGMACS weighted counts
disp("Fitting the four distributions to the BIGMACS weighted counts")
[mixLogNormPDF, ~, BM.MLN.fitInfo] = fitMixLogNorm(BM.weightedC, d.x, 2, d.S1.fitS.mlnReps, BM.numCpairs);
BM.MLN.nSR.x  = mixLogNormPDF(:,1);
BM.MLN.nSR.px = mixLogNormPDF(:,2);
[logNormPDF, ~, BM.LN.fitInfo] = fitMixLogNorm(BM.weightedC, d.x, 1, d.S1.fitS.mlnReps, BM.numCpairs);
BM.LN.nSR.x   = logNormPDF(:,1);
BM.LN.nSR.px  = logNormPDF(:,2);
[BM.Gam.nSR.px, BM.Gam.fitInfo] = fitGamma(BM.weightedC, d.x, BM.numCpairs);
BM.Gam.nSR.x  = d.x';
[BM.invGam.nSR.px, ~, BM.invGam.fitInfo] = fitInvGamma(BM.weightedC, d.x, BM.numCpairs);
BM.invGam.nSR.x = d.x';

% Transform the BIGMACS fits to log(nSR)
[BM.LN.lnSR.x, BM.LN.lnSR.px]         = px_to_pfx(BM.LN.nSR.x, BM.LN.nSR.px, logFH);
[BM.MLN.lnSR.x, BM.MLN.lnSR.px]       = px_to_pfx(BM.MLN.nSR.x, BM.MLN.nSR.px, logFH);
[BM.Gam.lnSR.x, BM.Gam.lnSR.px]       = px_to_pfx(BM.Gam.nSR.x, BM.Gam.nSR.px, logFH);
[BM.invGam.lnSR.x, BM.invGam.lnSR.px] = px_to_pfx(BM.invGam.nSR.x, BM.invGam.nSR.px, logFH);

fC = figure;
fC.Position = fig4Position;
panelLayout = tiledlayout(fig4Layout(1), fig4Layout(2), 'TileSpacing', 'tight', 'Padding', 'compact');
panelAxes   = gobjects(numel(fitMethodNames), 1);
for iMethod = 1:numel(fitMethodNames)
    % Get this method's fit struct (BM for BIGMACS, otherwise from the fit file)
    if iMethod == 1
        methodFit = BM;
    else
        methodFit = d.S1.(fitMethodNames(iMethod));
    end

    % Get the pooled histogram, scaled to cm per run, and the factor that
    % puts the pdfs on the same scale
    method_logHC = histcounts(log(methodFit.weightedC), logBinEdges);
    method_logHC = method_logHC.*methodFit.sedLength./numel(methodFit.weightedC);
    normF        = 1./sum(method_logHC.*bw);

    % Plot the histogram with the four fits and the published BIGMACS
    % lognormal. The fits are drawn G, IG, MLN, LN so that LN ends up on top;
    % the legend lists them in the order LN, MLN, G, IG from the handles
    panelAxes(iMethod) = nexttile;
    hold on
    box on
    histogram('BinCounts', method_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 1)
    gamLine    = plot(methodFit.Gam.lnSR.x, methodFit.Gam.lnSR.px.*(1/normF), 'LineWidth', 1.5, 'DisplayName', 'G', 'Color', fitColours.Gam, 'LineStyle', fitLineStyles.Gam);
    invGamLine = plot(methodFit.invGam.lnSR.x, methodFit.invGam.lnSR.px.*(1/normF), 'LineWidth', 1.5, 'DisplayName', 'IG', 'Color', fitColours.invGam, 'LineStyle', fitLineStyles.invGam);
    mlnLine    = plot(methodFit.MLN.lnSR.x, methodFit.MLN.lnSR.px.*(1/normF), 'LineWidth', 1.5, 'DisplayName', 'MLN', 'Color', fitColours.MLN, 'LineStyle', fitLineStyles.MLN);
    lnLine     = plot(methodFit.LN.lnSR.x, methodFit.LN.lnSR.px.*(1/normF), 'LineWidth', 1.5, 'DisplayName', 'LN', 'Color', fitColours.LN, 'LineStyle', fitLineStyles.LN);
    bmLine     = plot(BM.lnSR.x, BM.lnSR.px.*(1/normF), 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'BIGMACS');
    if iMethod == 1
        legendLines = [lnLine, mlnLine, gamLine, invGamLine, bmLine];
    end
    xlim(fig4Xlim)
    text(0.02, 0.9, panelLetters(iMethod) + ": " + panelNames(iMethod), 'Units', 'normalized', ...
        'FontSize', 12, 'FontWeight', 'bold')

    % Label the x axis on the bottom row and the y axis on the left column
    % only, and drop the tick labels elsewhere (all panels share their limits)
    panelRow = ceil(iMethod./fig4Layout(2));
    panelCol = mod(iMethod-1, fig4Layout(2)) + 1;
    if panelRow == fig4Layout(1)
        xlabel("log(NSR)")
    else
        xticklabels([])
    end
    if panelCol == 1
        ylabel("cm per run")
    else
        yticklabels([])
    end
end

% Put one shared legend along the bottom of the figure, in the same font
% size as the axes
sharedLegend             = legend(legendLines, 'Orientation', 'horizontal');
sharedLegend.Layout.Tile = 'south';
sharedLegend.FontSize    = 12;

% Make all the y limits the same, set by the tallest panel
maxYlim = 0;
for i = 1:numel(panelAxes)
    maxYlim = max(maxYlim, max(ylim(panelAxes(i))));
end
for i = 1:numel(panelAxes)
    ylim(panelAxes(i), [0 maxYlim.*1.05])
end

%% ===== Supplementary plots and tables =====
%% Individual-run histograms: median with 68% bounds, for log(nSR) and dt (was Figure 3)
% One row per method: BSamp, RSR0, RSR500, RSR1000 (individual-run fits)
irNames      = [bsampName, rsr0Name, rsr500Name, rsr1000Name] + "IR";
panelLetters = ["A", "B", "C", "D", "E", "F", "G", "H"];
numruns      = numel(d.S1.(irNames(1)).OneRunDatas);

% Set up the percentile rows (median and 68% bounds) of the sorted runs
medianRow = round(numruns*0.5);
lowRow    = round(numruns*0.16);
highRow   = round(numruns*0.84);

fIR = figure;
fIR.Position = [100, 100, 1000, 800]; % [left, bottom, width, height]
for iMethod = 1:numel(irNames)
    irFit = d.S1.(irNames(iMethod));

    % Plot the median log(nSR) histogram with 68% bounds
    subplot(4,2,2*iMethod-1)
    runHistCounts = sort(irFit.lnSRHistCounts, 1);
    hold on
    box on
    histogram('BinCounts', runHistCounts(medianRow, :), 'BinEdges', logBinEdges, 'FaceColor', [0.8 0.8 0.8])
    errorbar(logBinCenters, runHistCounts(medianRow, :), runHistCounts(medianRow, :)-runHistCounts(lowRow, :), ...
        runHistCounts(highRow, :)-runHistCounts(medianRow, :), 'LineStyle', 'none', 'DisplayName', '68% bounds')
    xlim([-2.5 2.5])
    xlabel("log(NSR)")
    ylabel("cm per run")
    text(0.02, 0.9, panelLetters(2*iMethod-1), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

    % Plot the median dt histogram with 68% bounds
    subplot(4,2,2*iMethod)
    runDtCounts = sort(irFit.agediffsWbC', 1);
    hold on
    box on
    histogram('BinCounts', runDtCounts(medianRow, :), 'BinEdges', dtBinEdges, 'FaceColor', [0.75 0.75 0.75])
    errorbar(dtBinEdges(1:end-1)+0.05, runDtCounts(medianRow, :), runDtCounts(medianRow, :)-runDtCounts(lowRow, :), ...
        runDtCounts(highRow, :)-runDtCounts(medianRow, :), 'LineStyle', 'none', 'DisplayName', '68% bounds')
    xlim([0 6])
    xlabel("\Deltat (kyr)")
    ylabel("cm per run")
    text(0.02, 0.9, panelLetters(2*iMethod), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
end

% Make all the y limits the same within each column of plots
% (results_Apr2026 set them from the last panel rather than the tallest;
% here the tallest panel sets them, as the comment there intended)
for iCol = 1:2
    maxYlim = 0;
    for i = iCol:2:8
        subplot(4,2,i)
        maxYlim = max(maxYlim, max(ylim));
    end
    for i = iCol:2:8
        subplot(4,2,i)
        ylim([0 maxYlim.*1.05])
    end
end

%% Table S1: preferred distribution (lowest BICeff) across the individual-run fits
% Table 2 ranks the four distributions by BICeff for the pooled fits; this
% does the same for every individual run (one age-model realisation per
% core) of each IR method and gives the percentage of runs in which each
% distribution had the lowest BICeff. Runs where any of the four fits failed
% (NaN BICeff) are counted in numRunsFailed and left out of the percentages
irPreferredPct = NaN(numel(irNames), numel(distNames));
numRunsFitted  = NaN(numel(irNames), 1);
numRunsFailed  = NaN(numel(irNames), 1);
for iMethod = 1:numel(irNames)
    % Get the BICeff of each distribution in each run (runs x dists)
    IR         = d.S1.(irNames(iMethod));
    numIRruns  = numel(IR.LN.fits);
    irBICeff   = NaN(numIRruns, numel(distNames));
    for iDist = 1:numel(distNames)
        for i = 1:numIRruns
            irBICeff(i, iDist) = IR.(distNames(iDist)).fits{i}.BICtaeheefix;
        end
    end

    % Leave out runs where any fit failed, so every counted run compared all four
    failedLog              = any(isnan(irBICeff), 2);
    irBICeff               = irBICeff(~failedLog, :);
    numRunsFitted(iMethod) = sum(~failedLog);
    numRunsFailed(iMethod) = sum(failedLog);

    % Find the lowest BICeff in each run and count how often each distribution won
    [~, preferredInd] = min(irBICeff, [], 2);
    for iDist = 1:numel(distNames)
        irPreferredPct(iMethod, iDist) = 100.*sum(preferredInd == iDist)./numRunsFitted(iMethod);
    end
end
IRpreferred_table = array2table(round(irPreferredPct, 1), 'VariableNames', distNames, 'RowNames', irNames);
IRpreferred_table.numRunsFitted = numRunsFitted;
IRpreferred_table.numRunsFailed = numRunsFailed

%% Pooled RSR500 histogram with the LN and IG fits on an inverse nSR axis
% The RSR500 pooled (AR) and individual-run (IR) fits used by this and the
% next supplementary section
AR = d.S1.(rsr500Name);
IR = d.S1.(rsr500Name + "IR");

invBinEdges = 0:0.1:6;
f_inv       = @(x)(1./x);
[invLN.x, invLN.px]         = px_to_pfx(AR.LN.nSR.x, AR.LN.nSR.px, f_inv);
[invInvGam.x, invInvGam.px] = px_to_pfx(AR.invGam.nSR.x, AR.invGam.nSR.px, f_inv);

% Get the pooled inverse nSR histogram, scaled to cm per run
AR_invHC = histcounts(1./AR.weightedC, invBinEdges);
AR_invHC = AR_invHC.*AR.sedLength./numel(AR.weightedC);
normF    = 1./sum(AR_invHC.*0.1);

figure;
hold on
histogram('BinCounts', AR_invHC, 'BinEdges', invBinEdges, 'HandleVisibility', 'off', 'FaceColor', [0.75 0.75 0.75])
plot(invLN.x, invLN.px*(1/normF), 'LineWidth', 1, 'DisplayName', 'LN', 'Color', fitColours.LN, 'LineStyle', fitLineStyles.LN)
plot(invInvGam.x, invInvGam.px*(1/normF), 'LineWidth', 1, 'DisplayName', 'IG', 'Color', fitColours.invGam, 'LineStyle', fitLineStyles.invGam)
xlim([0 6])
xlabel("inverse NSR")
ylabel("cm per run")
legend()

%% Fits to the pooled data vs fits to each individual run (RSR500)
% Evidence that fitting to the pooled data of many MC runs gives results
% representative of fitting to the data of each MC run
lbe       = -4:0.1:4;   % Log bin edges
numIRfits = numel(IR.LN.fits);

% Plot pooled histograms with the fit to each MC run and the fit to the pooled data
figure;
for iDist = 1:numel(distNames)
    subplot(2,4,iDist)
    ARvsIRplots_wph(AR.(distNames(iDist)), IR.(distNames(iDist)), AR, lbe, fitColours.(distNames(iDist)))
    text(0.02, 0.9, panelLetters(iDist), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
end

% Plot the fitted parameters of each MC run against the pooled fit: LN
subplot(2,4,5)
hold on
for i = 1:numIRfits
    plot(IR.LN.fits{i}.mu, sqrt(IR.LN.fits{i}.Sigma), '.k')
end
plot(AR.LN.fitInfo.mu, sqrt(AR.LN.fitInfo.Sigma), 'o', 'Color', fitColours.LN, 'MarkerFaceColor', fitColours.LN)
xlabel("Mu")
ylabel("Sigma")
text(0.02, 0.9, 'E', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% MLN (mean and standard deviation of the nSR pdf)
subplot(2,4,6)
hold on
for i = 1:numIRfits
    plot(IR.MLN.nSR.mu(i), sqrt(IR.MLN.nSR.var(i)), '.k')
end
plot(AR.MLN.nSR.mu, sqrt(AR.MLN.nSR.var), 'o', 'Color', fitColours.MLN, 'MarkerFaceColor', fitColours.MLN)
xlabel("Mean")
ylabel("Standard Deviation")
text(0.02, 0.9, 'F', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% Gam
subplot(2,4,7)
hold on
for i = 1:numIRfits
    plot(IR.Gam.fits{i}.alpha, IR.Gam.fits{i}.beta, '.k')
end
plot(AR.Gam.fitInfo.alpha, AR.Gam.fitInfo.beta, 'o', 'Color', fitColours.Gam, 'MarkerFaceColor', fitColours.Gam)
xlabel("Alpha")
ylabel("Beta")
text(0.02, 0.9, 'G', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% invGam
subplot(2,4,8)
hold on
for i = 1:numIRfits
    plot(IR.invGam.fits{i}.alpha, IR.invGam.fits{i}.beta, '.k')
end
plot(AR.invGam.fitInfo.alpha, AR.invGam.fitInfo.beta, 'o', 'Color', fitColours.invGam, 'MarkerFaceColor', fitColours.invGam)
xlabel("Alpha")
ylabel("Beta")
text(0.02, 0.9, 'H', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')

% Make the y limits of the histogram panels the same (tallest panel)
maxYlim = 0;
for i = 1:4
    subplot(2,4,i)
    maxYlim = max(maxYlim, max(ylim));
end
for i = 1:4
    subplot(2,4,i)
    ylim([0 3200])
end

%% Compare RSR500 distributions for the different weightings
% Load each weighting's fit (the main fit is already loaded)
weightFits = cell(numel(runS.weightings), 1);
weightFits{1} = d;
for iWeight = 2:numel(runS.weightings)
    weightFit = load(runS.fitFiles(iWeight), "d");
    weightFits{iWeight} = weightFit.d;
end

%Y axis labels for the  3 weightings
weight_yls = [];

% Plot the pooled RSR500 histogram and the four fits for each weighting,
% one panel per weighting, with one shared legend below the x axis label
figure;
weightLayout = tiledlayout(numel(runS.weightings), 1, 'TileSpacing', 'tight', 'Padding', 'compact');
weightAxes   = gobjects(numel(runS.weightings), 1);
for iWeight = 1:numel(runS.weightings)
    weightAR = weightFits{iWeight}.S1.(rsr500Name);
    weightAxes(iWeight) = nexttile;

    % Set the per-run total the histogram is scaled to, which depends on
    % the weighting: cm of sediment (depth), kyr of time (age, sedTimeSpan
    % is in yr) or number of age pairs (none)
    if runS.weightings(iWeight) == "depth"
        weightTotal = weightAR.sedLength;
        weightUnits = "cm per run";
    elseif runS.weightings(iWeight) == "age"
        weightTotal = weightAR.sedTimeSpan./1000;
        weightUnits = "kyr per run";
    else
        weightTotal = weightAR.numCpairs;
        weightUnits = "age pairs per run";
    end

    % Get the pooled histogram, scaled to that total per run (weightedC is
    % a sample of points, so its own count is arbitrary)
    weight_logHC = histcounts(log(weightAR.weightedC), logBinEdges);
    weight_logHC = weight_logHC.*weightTotal./numel(weightAR.weightedC);
    normF        = 1./sum(weight_logHC.*bw);
    histogram('BinCounts', weight_logHC, 'BinEdges', logBinEdges, 'FaceColor', [0.75 0.75 0.75], 'FaceAlpha', 1, 'HandleVisibility', 'off')
    hold on
    plot(weightAR.LN.lnSR.x, weightAR.LN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName', 'LN', 'Color', fitColours.LN, 'LineStyle', fitLineStyles.LN)
    plot(weightAR.MLN.lnSR.x, weightAR.MLN.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName', 'MLN', 'Color', fitColours.MLN, 'LineStyle', fitLineStyles.MLN)
    plot(weightAR.Gam.lnSR.x, weightAR.Gam.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName', 'G', 'Color', fitColours.Gam, 'LineStyle', fitLineStyles.Gam)
    plot(weightAR.invGam.lnSR.x, weightAR.invGam.lnSR.px*(1/normF), 'LineWidth', 1, 'DisplayName', 'IG', 'Color', fitColours.invGam, 'LineStyle', fitLineStyles.invGam)
    text(0.02, 0.85, panelLetters(iWeight), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
    ylabel(weightUnits)
    xlim([-2.5 2.5])

    % Keep the x tick labels on the bottom panel only
    if iWeight < numel(runS.weightings)
        xticklabels([])
    end
end
xlabel("log(NSR)")
weightLegend             = legend(weightAxes(1), 'Orientation', 'horizontal');
weightLegend.Layout.Tile = 'south';

% Make the BIC table for RSR500 under each weighting
bicWeightValues = NaN(numel(runS.weightings), numel(distNames));
for iWeight = 1:numel(runS.weightings)
    for j = 1:numel(distNames)
        bicWeightValues(iWeight,j) = weightFits{iWeight}.S1.(rsr500Name).(distNames(j)).fitInfo.BICtaeheefix;
    end
end
BIC_table_RSR500 = array2table(round(bicWeightValues), 'VariableNames', distNames, 'RowNames', runS.weightings)


function ARvsIRplots_wph(ARd, IRd, AR, logBinEdges, fitColour)
% Plot the pooled (AR) fit, in fitColour, on top of the individual-run (IR)
% fits, with the pooled histogram (cm per run) underneath
AR_logHC = histcounts(log(AR.weightedC), logBinEdges);
AR_logHC = AR_logHC.*AR.sedLength./numel(AR.weightedC);
normF    = 1./sum(AR_logHC.*uniquetol(diff(logBinEdges), 1e-5));
hold on
histogram('BinCounts', AR_logHC, 'BinEdges', logBinEdges, 'HandleVisibility', 'off', 'FaceColor', [0.8 0.8 0.8])
plot(IRd.lnSR.x, IRd.lnSR.px(:,1).*(1/normF), 'Color', [0 0 0 0.1], 'DisplayName', 'Individual Run')
plot(IRd.lnSR.x, IRd.lnSR.px(:,2:end).*(1/normF), 'Color', [0 0 0 0.1], 'HandleVisibility', 'off')
plot(ARd.lnSR.x, ARd.lnSR.px.*(1/normF), 'Color', fitColour, 'LineWidth', 1, 'DisplayName', 'Pooled')
xlim([-2.5 2.5])
xlabel("log(NSR)")
ylabel("cm per run")
end
