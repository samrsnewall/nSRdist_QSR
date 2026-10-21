%% samplingDirectionSensitivity_review1_allcores - evidence 2 for the response to Reviewer 1
%
% Compares the RSR500 histograms and fitted distributions of the three
% sequential sampling orders (shallow2deep, deep2shallow and random) on the
% whole production core set. The calcData file for runTag already exists;
% this script fits it with fitData on every core (commonCoresMethods empty),
% saving the fit file in fitFolder so it doesn't collide with other scripts'
% fits of the same calcData file. With plotsOnly = true the fit file in
% fitFolder is loaded and plotted without refitting.
%
% fitData is run as a script in this workspace and picks up the fitOverride
% struct set below; every other setting (the method lists, the fit
% options) comes from the settings block inside fitData. The calcData
% section is kept, commented out, as a record of how the file was made.
%
% The script cds to the repository root because calcData and fitData use
% paths relative to it.

cd(fileparts(fileparts(mfilename('fullpath'))))

%% Settings
runTag            = "QSR";              % Used in the calcData and fit file names
coreSelection     = "good";   % S.coreSelection: "good" for the full production core set, or a string array of core names
maxNumCores       = [];       % S.maxNumCores: keep only the first this-many chosen cores alphabetically ([] = all)
numruns           = 10000;    % S.numruns: Monte Carlo iterations for BSamp, RSRx and rejection
mainSeed          = 1;        % rng seed before the main calcData run
noiseSeed         = 2;        % rng seed before the rejection-only run (the noise baseline)
weighting         = "depth";  % fitS.weighting: "depth", "age" or "none"; the fit file suffix follows from it
maxFitPoints      = 5e6;      % fitS.maxFitPoints: cap on weighted points passed to the fitters. With all cores at
                              %   10000 runs the weighted data is ~1.6e8 points and each pooled fit takes hours;
                              %   5e6 gives the same fits in minutes
weightingPipeline = "sample"; % fitS.weightingPipeline: "replicate" or "sample" (same distribution, much faster)
commonCoresMethods = []; % fitS.commonCoresMethods: every fit uses only cores all of these produced results for ([] = all cores)
coreFraction       = 1;         % fitS.coreFraction: random fraction of those cores to fit (1 = all, 0.5 = a random half)
coreFractionSeed   = 5;         % fitS.coreFractionSeed 
fitFolder         = "Results/results_QSR"; % Where this script's fit file is saved (fitOverride.saveFolder).
                              %   The manuscript fit file is read directly, because its fit settings (every
                              %   core, maxFitPoints, weightingPipeline, weighting) are identical to the ones here
plotsOnly         = true;    % If true, fitData is not run: the fit file already in fitFolder is loaded and plotted

% Set up the file names; fitSuffix must match the rule in fitData.m
mainFile    = "Results/dataT_" + runTag;
fitSuffix   = "_fit_" + weighting + "weight";
mainFitFile = fullfile(fitFolder, "dataT_" + runTag + fitSuffix + ".mat");

% calcData and fitData define many short variable names, so the workspace is
% cleared of everything except these between runs
keepVars = {'runTag' 'coreSelection' 'maxNumCores' 'numruns' 'mainSeed' 'noiseSeed' 'weighting' ...
            'maxFitPoints' 'weightingPipeline' 'fitFolder' 'plotsOnly' 'mainFile' 'fitSuffix' 'mainFitFile' ...
            'commonCoresMethods' 'coreFraction' 'coreFractionSeed' 'keepVars'};

%% Check the fit file: it must exist to plot only, and must not exist to fit
if plotsOnly && ~isfile(mainFitFile)
    error("Fit file not found: " + mainFitFile + " - set plotsOnly = false to make it")
elseif ~plotsOnly && isfile(mainFitFile)
    error("Fit file already exists: " + mainFitFile + " - change fitFolder or delete the file, or set plotsOnly = true to plot it")
end
if ~isfolder(fitFolder)
    mkdir(fitFolder)
end

%% Run calcData: every sampling order and the rejection sampler
% if ~fitsOnly
%     clearvars('-except', keepVars{:})
%     calcOverride.dataOutputFile = "dataT_" + runTag;
%     calcOverride.coreSelection  = coreSelection;
%     calcOverride.maxNumCores    = maxNumCores;
%     calcOverride.numruns        = numruns;
% 
%     rng(mainSeed)
%     disp("===== calcData (main run) started " + string(datetime("now")))
%     calcData
%     disp("===== calcData (main run) finished " + string(datetime("now")))
% end
% 

%% Fit the main run
if ~plotsOnly
    clearvars('-except', keepVars{:})
    fitOverride.filepath           = mainFile;
    fitOverride.saveFolder         = fitFolder;
    fitOverride.weighting          = weighting;
    fitOverride.maxFitPoints       = maxFitPoints;
    fitOverride.weightingPipeline  = weightingPipeline;
    fitOverride.commonCoresMethods = commonCoresMethods;
    fitOverride.coreFraction       = coreFraction;
    fitOverride.coreFractionSeed   = coreFractionSeed;

    % Fit only the RSR columns (fitData always fits every RSR column): no
    % Bchron or rejection fits, since this script only tests the orders
    fitOverride.lin2014Methods     = strings(1, 0);
    fitOverride.unthinnedMethods   = strings(1, 0);

    disp("===== fitData (main run) started " + string(datetime("now")))
    fitData
    disp("===== fitData (main run) finished " + string(datetime("now")))
    clearvars('-except', keepVars{:})
end

%% Load the fit results
mainFit  = load(mainFitFile, "d");
dMain    = mainFit.d;

%% Histogram and figure settings
bw             = 0.1;              % Bin width in log(nSR) space
logBinEdges    = -5:bw:5;
logBinCenters  = (logBinEdges(1:end-1) + logBinEdges(2:end)).*0.5;
xlimLog        = [-2.5 2.5];       % x-limits for all log(NSR) panels
chosenOrder    = "";               % Tag of the order used in the manuscript (shallow2deep, stored with no suffix); residuals are taken relative to it
rsrMinAgeDiffs = [0, 500, 1000];   % The restrictions calcData runs for each sampling order
distNames      = ["LN", "MLN", "Gam", "invGam"];          % Field names in the fit structs
distLabels     = ["Lognormal", "Mixture Lognormal", "Gamma", "Inverse Gamma"];
panelLetters   = ["A", "B", "C", "D", "E", "F", "G", "H"];
orderLines     = ["-", "--", ":"];  % Line styles for the chosen order, then the alternatives

% Set up the sampling-order tags used in the fit names, chosen order first
allOrderTags   = replace(dMain.S.sampleOrders, ["shallow2deep", "deep2shallow"], ["", "d2s"]);
otherOrderTags = allOrderTags(allOrderTags ~= chosenOrder);
orderTags      = [chosenOrder, otherOrderTags];

% Get the weight inflator used when the fitted data were replicated
weightInflator = dMain.S1.fitS.PooledRuns.weightRepInflator;

%% Calculate the log(nSR) histogram of every fit
% HC holds the counts; P holds them normalised by their own total, so the
% differences below show shape only, not differences in the number of counts
fitNames = string(fieldnames(dMain.S1));
fitLog   = ~ismember(fitNames, ["fitS", "chooseLog", "numCores", "chosenCoreNames"]) & ~endsWith(fitNames, "IR");
fitNames = fitNames(fitLog);
for iFit = 1:numel(fitNames)
    fitNSRs    = dMain.S1.(fitNames(iFit)).weightedC;
    histCounts = histcounts(log(fitNSRs), logBinEdges)./weightInflator;
    HC.(fitNames(iFit)) = histCounts;
    P.(fitNames(iFit))  = histCounts./sum(histCounts);
end

%% ===== Reviewer figures =====
%% Fig 1: RSR500 log(NSR) histograms, one panel per sampling order
figure;
for k = 1:numel(orderTags)
    subplot(numel(orderTags), 1, k)
    histogram('BinCounts', HC.("RSR500" + orderTags(k)), 'BinEdges', logBinEdges, 'FaceColor', [0.8 0.8 0.8]);
    text(0.02, 0.88, panelLetters(k), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
    text(0.98, 0.88, "RSR500" + orderTags(k), 'Units', 'normalized', 'HorizontalAlignment', 'right')
    box on
    xlim(xlimLog)
    ylabel("cm")
end
xlabel("log(NSR)")

% Set the same y-axis limits on all panels
ax = findobj(gcf, 'Type', 'axes');
yl = cell2mat(get(ax, 'YLim'));
set(ax, 'YLim', [0 max(yl(:,2))]);

%% Fig 2: residuals of the alternative sampling orders relative to the chosen one
figure;
for k = 1:numel(otherOrderTags)
    subplot(numel(otherOrderTags), 1, k)
    hold on
    histDiff = P.("RSR500" + otherOrderTags(k)) - P.("RSR500" + chosenOrder);
    bar(logBinCenters, histDiff, 1, 'FaceColor', [0.8 0.8 0.8]);
    plot(xlimLog, [0 0], 'Color', 'k', 'LineWidth', 0.5)
    text(0.02, 0.88, panelLetters(k), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
    text(0.98, 0.88, "RSR500" + otherOrderTags(k) + " - RSR500" + chosenOrder, 'Units', 'normalized', 'HorizontalAlignment', 'right')
    box on
    xlim(xlimLog)
    ylabel("fraction of counts")
end
xlabel("log(NSR)")

% Set the same symmetric y-axis limits on all panels
ax = findobj(gcf, 'Type', 'axes');
yl = cell2mat(get(ax, 'YLim'));
ymax = max(abs(yl(:)));
set(ax, 'YLim', [-ymax ymax]);

%% Fig 3: fitted log(NSR) pdfs, one panel per distribution, one line per sampling order
figure;
for iDist = 1:numel(distNames)
    subplot(2, 2, iDist)
    hold on
    for k = 1:numel(orderTags)
        fitPdf = dMain.S1.("RSR500" + orderTags(k)).(distNames(iDist)).lnSR;
        plot(fitPdf.x, fitPdf.px, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', orderLines(k), 'DisplayName', "RSR500" + orderTags(k))
    end
    text(0.02, 0.92, panelLetters(iDist), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
    title(distLabels(iDist), 'FontWeight', 'normal')
    box on
    xlim(xlimLog)
    xlabel("log(NSR)")
    ylabel("pdf")
    if iDist == 1
        legend('Location', 'northwest')
    end
end

% Set the same y-axis limits on all panels
ax = findobj(gcf, 'Type', 'axes');
yl = cell2mat(get(ax, 'YLim'));
set(ax, 'YLim', [0 max(yl(:,2))]);

%% Fig 4: residuals of the fitted pdfs relative to the chosen sampling order
figure;
for iDist = 1:numel(distNames)
    subplot(2, 2, iDist)
    hold on
    chosenPdf = dMain.S1.("RSR500" + chosenOrder).(distNames(iDist)).lnSR;
    for k = 1:numel(otherOrderTags)
        fitPdf = dMain.S1.("RSR500" + otherOrderTags(k)).(distNames(iDist)).lnSR;
        plot(fitPdf.x, fitPdf.px - chosenPdf.px, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', orderLines(k+1), 'DisplayName', "RSR500" + otherOrderTags(k) + " - RSR500" + chosenOrder)
    end
    plot(xlimLog, [0 0], 'Color', 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off')
    text(0.02, 0.92, panelLetters(iDist), 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold')
    title(distLabels(iDist), 'FontWeight', 'normal')
    box on
    xlim(xlimLog)
    xlabel("log(NSR)")
    ylabel("pdf difference")
    if iDist == 1
        legend('Location', 'northwest')
    end
end

% Set the same symmetric y-axis limits on all panels
ax = findobj(gcf, 'Type', 'axes');
yl = cell2mat(get(ax, 'YLim'));
ymax = max(abs(yl(:)));
set(ax, 'YLim', [-ymax ymax]);


%% ===== Reviewer figures: RSRx sampling orders on every core =====
% Two figures that summarise the order test for the response to Reviewer 1.
% Histograms are shown as probability density in log(nSR) (fraction of
% weighted counts per bin divided by the bin width), so that the histogram
% and its residuals are in the same units as the fitted pdfs in the second
% figure. Only the three sequential orders are shown: the joint (rejection)
% sampler gave up on some of these cores, so it is compared on the cores it
% could handle in samplingDirectionSensitivity_review1_jointSampleCores
reviewerMinAgeDiff = 500;                                        % Minimum age difference of the RSRx fits to show (0, 500 or 1000 yr)
reviewerNames   = "RSR" + reviewerMinAgeDiff + [chosenOrder, "d2s", "random"];
reviewerLabels  = ["shallow-to-deep (manuscript)", "deep-to-shallow", "random order"];
reviewerColors  = ['k', 'r', 'b'];
reviewerLines   = ["-", "--", ":"];
reviewerWidths  = [1.8, 1.2, 1.4];
reviewerYmax    = 0.9;                                           % y-limit for the density histogram panel
reviewerFitYmax = 1;                                             % y-limit for the fitted pdf panels (room for the legend)
reviewerResid   = 0.045;                                         % Symmetric y-limit for the residual panel
reviewerFigWidth   = 16;                                         % Figure width (cm), to fit the text width of the response letter
reviewerFig5Height = 11;                                         % Height of Fig 5 (cm)
reviewerFig6Height = 13;                                         % Height of Fig 6 (cm)
reviewerFontSize   = 10;                                         % Font size (pt) of tick labels and legend; axis labels are 1.1x this
reviewerLetterSize = 12;                                         % Font size (pt) of the bold panel letters

% Convert the histograms to density
for iName = 1:numel(reviewerNames)
    dens.(reviewerNames(iName)) = P.(reviewerNames(iName))./bw;
end

%% Fig 5: RSRx density histograms, with the residuals from the chosen order underneath
% Uses tiledlayout rather than subplot so that the one legend can sit in a
% row along the bottom of the figure ('south' tile)
figure('Units', 'centimeters', 'Position', [2 2 reviewerFigWidth reviewerFig5Height]);
tiledlayout(3, 1, 'TileSpacing', 'none', 'Padding', 'compact')
densAxes = nexttile([2 1]);
hold on
for iName = 1:numel(reviewerNames)
    stairs(logBinEdges(1:end-1), dens.(reviewerNames(iName)), 'Color', reviewerColors(iName), ...
        'LineStyle', reviewerLines(iName), 'LineWidth', reviewerWidths(iName), 'DisplayName', reviewerLabels(iName))
end
text(0.02, 0.92, panelLetters(1), 'Units', 'normalized', 'FontSize', reviewerLetterSize, 'FontWeight', 'bold')
box on
xlim(xlimLog)
ylim([0 reviewerYmax])
set(gca, 'XTickLabel', [])
ylabel("probability density")

% Plot the residuals relative to the chosen order
nexttile
hold on
for iName = 2:numel(reviewerNames)
    densDiff = dens.(reviewerNames(iName)) - dens.(reviewerNames(1));
    stairs(logBinEdges(1:end-1), densDiff, 'Color', reviewerColors(iName), ...
        'LineStyle', reviewerLines(iName), 'LineWidth', reviewerWidths(iName), 'HandleVisibility', 'off')
end
plot(xlimLog, [0 0], 'Color', 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off')
text(0.02, 0.85, panelLetters(2), 'Units', 'normalized', 'FontSize', reviewerLetterSize, 'FontWeight', 'bold')
box on
xlim(xlimLog)
ylim([-reviewerResid reviewerResid])
yticks(-0.02:0.02:0.02)                                          % No tick at the top edge, so it doesn't clash with the 0 above
xlabel("log(NSR)")
ylabel(["density difference", "from shallow-to-deep"])

% Set the font size on both panels (fixing it also stops tiledlayout shrinking it)
set(findobj(gcf, 'Type', 'axes'), 'FontSize', reviewerFontSize)

% Put one horizontal legend along the bottom of the figure, from the series of the density panel
densLegend = legend(densAxes, 'Orientation', 'horizontal', 'Box', 'on', 'FontSize', reviewerFontSize);
densLegend.Layout.Tile = 'south';

%% Fig 6: fitted log(NSR) pdfs of the four distributions, one panel each
figure('Units', 'centimeters', 'Position', [2 2 reviewerFigWidth reviewerFig6Height]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
for iDist = 1:numel(distNames)
    nexttile
    hold on
    for iName = 1:numel(reviewerNames)
        fitPdf = dMain.S1.(reviewerNames(iName)).(distNames(iDist)).lnSR;
        plot(fitPdf.x, fitPdf.px, 'Color', reviewerColors(iName), 'LineStyle', reviewerLines(iName), ...
            'LineWidth', reviewerWidths(iName), 'DisplayName', reviewerLabels(iName))
    end
    % Label the panel: bold letter top left, distribution name top right. Two
    % separate texts, because switching \bf to \rm in one string widens the spaces
    text(0.02, 0.92, panelLetters(iDist), 'Units', 'normalized', 'FontSize', reviewerLetterSize, 'FontWeight', 'bold')
    text(0.98, 0.92, distLabels(iDist), 'Units', 'normalized', 'FontSize', reviewerFontSize, 'HorizontalAlignment', 'right')
    box on
    xlim(xlimLog)
    ylim([0 reviewerFitYmax])
    xlabel("log(NSR)")
    ylabel("fitted pdf")
end

% Set the font size on all four panels (fixing it also stops tiledlayout shrinking it)
set(findobj(gcf, 'Type', 'axes'), 'FontSize', reviewerFontSize)

% Put one horizontal legend along the bottom of the figure
fitLegend = legend('Orientation', 'horizontal', 'Box', 'on', 'FontSize', reviewerFontSize);
fitLegend.Layout.Tile = 'south';

%% ===== BICeff ranking check =====
% The manuscript's family choice rests on the BICeff ranking of the four
% distributions. Print the ranking of every fit (lowest BICeff first) with
% the gap between the best and the second best, then check that each
% alternative sampling order ranks them the same way as the chosen order

% Print the ranking of every fit in the main run
disp("BICeff ranking (best first), and the gap from 1st to 2nd:")
for iFit = 1:numel(fitNames)
    fitName = fitNames(iFit);
    BICeff  = NaN(1, numel(distNames));
    for iDist = 1:numel(distNames)
        BICeff(iDist) = dMain.S1.(fitName).(distNames(iDist)).fitInfo.BICtaeheefix;
    end
    [sortedBICeff, rankOrder] = sort(BICeff);
    rankOrders.(fitName) = rankOrder;

    rankedText = strjoin(distNames(rankOrder) + " (" + compose("%.1f", sortedBICeff) + ")", " > ");
    bicGap     = sortedBICeff(2) - sortedBICeff(1);
    disp("   " + fitName + ": " + rankedText + "   gap " + compose("%.1f", bicGap))
end

% Check each alternative sampling order against the chosen one, for every restriction
for iDt = 1:numel(rsrMinAgeDiffs)
    chosenName = "RSR" + rsrMinAgeDiffs(iDt) + chosenOrder;
    for k = 1:numel(otherOrderTags)
        otherName = "RSR" + rsrMinAgeDiffs(iDt) + otherOrderTags(k);
        if isequal(rankOrders.(otherName), rankOrders.(chosenName))
            disp(otherName + " ranks the distributions the same as " + chosenName)
        else
            warning(otherName + " ranks the distributions differently from " + chosenName)
        end
    end
end
