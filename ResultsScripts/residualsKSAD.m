%% residualsKSAD - residual histograms and KS / AD distances (Reviewer 1, point 2)
%
% Two more of the tail-focused GOF checks the reviewer asked for, run on the
% same fit file as qq_tailGOF and quantileSegmentFits.
%e
% Residual histograms: bar plots of the per-bin difference between the
% weighted empirical distribution and each fitted family, in ln(nSR) bins.
% The bars are standardised (Pearson) residuals, (observed - expected) /
% standard error of the expected fraction at numEff pairs, so every bin is
% on a "how many standard errors off" scale and tail misfit is as visible
% as body misfit. Raw probability differences (as in
% samplingDirectionSensitivity) are available with resType = "raw"; they
% are near zero in the tails by construction, because both the observed and
% expected fractions are small there.
%
% KS and AD: the Kolmogorov-Smirnov distance D (largest gap between the
% weighted empirical CDF and the fitted CDF) and the Anderson-Darling
% statistic A2, which weights the squared CDF gap by 1/(F(1-F)) and so
% emphasises the tails. Both are reported as descriptive distances, not
% with textbook p-values: those would assume independent observations, the
% same assumption Reviewer 1 objected to for BICeff, and would also ignore
% that the parameters were fitted to the same data. D is also given within
% each quantile segment, which is the "segment-wise error" the reviewer
% mentioned. numEff = numCpairs is used wherever a sample size enters, as
% for BICeff.
%
% Run from the ResultsScripts folder.

clear; close all
addpath('../Functions')

%% Settings
fitFile   = "../Results/results_QSR/dataT_QSR_fit_depthweight.mat";
families  = ["LN", "MLN", "Gam", "invGam"];
familyLbl = ["Lognormal", "Mixture Lognormal", "Gamma", "Inverse Gamma"];
segPct    = [0.05 0.95];    % Empirical percentiles that separate the three quantile segments
resType   = "raw";      % Residual scale for the bar plots: "pearson" (standardised at numEff) or "raw" (probability difference)
resBW     = 0.2;            % Residual histogram bin width in ln(nSR)
minExpCount = 0;            % Bins whose expected count at numEff is below this are drawn hollow (too sparse to read)
seLines   = 2;              % Reference lines at +-seLines standard errors on the Pearson plots
pearsonYLim = 10;           % y-axis limits of +-pearsonYLim standard errors on the Pearson plots; bars beyond
                            % this are clipped and marked with a triangle at the axis edge (invGam is far off)
saveFig   = true;
saveTable = true;
dateTag   = string(datetime("today", "Format", "ddMMMyy"));

% Rows of the residual figure, top to bottom (manuscript order). RSR0 has no
% minimum dt and a much heavier upper tail, so its row uses the wider ln(nSR)
% range; the other rows share the narrow one
methods       = ["BMedianLin", "BSampLin", "RSR0", "RSR500", "RSR1000"];
methodLbl     = ["BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];
wideMethods   = "RSR0";
resRange      = [-3 3];     % Residual histogram range in ln(nSR) for most rows
resRangeWide  = [-3 6];     % ... and for wideMethods
figName       = "../Figures/Residuals_" + dateTag + ".png";
size_cm       = [26 27];

tableFile = "../Results/KS_AD_distances_" + dateTag + ".csv";

% Paper tables: one for KS D and one for AD A2, methods as rows and families as columns.
% paperRows sets the row order (labels as in methodLbl); paperFamLbl gives the
% column headings, as the families are named in the manuscript
paperRows   = ["BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];
paperFamLbl = ["Lognormal", "Mixture Lognormal", "Gamma", "Inverse Gamma"];
KStableFile = "../Results/KS_table_" + dateTag + ".csv";
ADtableFile = "../Results/AD_table_" + dateTag + ".csv";

%% Load the fit results
load(fitFile, "d")

%% Residual histograms
numMethods = numel(methods);
numFams    = numel(families);

figure('Units', 'centimeters', 'Position', [2 2 size_cm]);
tl = tiledlayout(numMethods, numFams, 'TileSpacing', 'compact', 'Padding', 'compact');
resAxes = gobjects(numMethods, numFams);

for i = 1:numMethods
    % Set the residual histogram bins for this row
    if ismember(methods(i), wideMethods)
        resEdges = resRangeWide(1):resBW:resRangeWide(2);
    else
        resEdges = resRange(1):resBW:resRange(2);
    end
    resCenters = (resEdges(1:end-1) + resEdges(2:end)).*0.5;

    % Get the weighted data for this method and its effective sample size
    methodFit = d.S1.(methods(i));
    weightedC = methodFit.weightedC;
    if iscell(weightedC), weightedC = cell2mat(weightedC(:)); end
    weightedC = weightedC(isfinite(weightedC) & weightedC > 0);
    lnSRvals  = log(weightedC);
    numEff    = round(mean(methodFit.numCpairs));

    % Calculate the observed fraction of pairs in each ln(nSR) bin
    binCounts = histcounts(lnSRvals, resEdges);
    pObs      = binCounts./numel(lnSRvals);

    for j = 1:numFams
        % Calculate the fraction each fitted family puts in the same bins
        cdfAtEdges = fittedCDF(methodFit.(families(j)).fitInfo, families(j), exp(resEdges));
        pExp       = diff(cdfAtEdges);

        % Calculate the residual on the chosen scale
        if resType == "pearson"
            seExp     = sqrt(pExp.*(1 - pExp)./numEff);
            residuals = (pObs - pExp)./seExp;
        else
            residuals = pObs - pExp;
        end

        % Find bins too sparse to read (expected count below minExpCount)
        sparseLog = numEff.*pExp < minExpCount;

        % Split the residuals into filled (well-populated) and hollow (sparse)
        % bars; NaNs keep the bar spacing right when bins are skipped
        residualsFilled = residuals;
        residualsHollow = residuals;
        residualsFilled(sparseLog)  = NaN;
        residualsHollow(~sparseLog) = NaN;

        % Plot the residual bars
        resAxes(i,j) = nexttile; hold on; box on
        bar(resCenters, residualsFilled, 1, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
        bar(resCenters, residualsHollow, 1, 'FaceColor', 'none', 'EdgeColor', [0.6 0.6 0.6]);
        yline(0, 'k-');
        if resType == "pearson"
            yline( seLines, ':', 'Color', [0.3 0.3 0.3]);
            yline(-seLines, ':', 'Color', [0.3 0.3 0.3]);

            % Mark bars that go beyond the fixed y-limits with a triangle at the edge
            clippedHighLog = residuals >  pearsonYLim;
            clippedLowLog  = residuals < -pearsonYLim;
            plot(resCenters(clippedHighLog),  0.95.*pearsonYLim.*ones(1, sum(clippedHighLog)), '^', 'MarkerSize', 4, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none')
            plot(resCenters(clippedLowLog),  -0.95.*pearsonYLim.*ones(1, sum(clippedLowLog)),  'v', 'MarkerSize', 4, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none')
        end
        xlim([-3 3])

        % Label the grid
        if i == 1, title(familyLbl(j), 'FontWeight', 'normal'); end
        if j == 1
            if resType == "pearson"
                ylabel(methodLbl(i) + newline + "(obs - exp) / SE")
            else
                ylabel(methodLbl(i))
            end
        end
        if i == numMethods, xlabel("ln(nSR)"); end
    end
end

% Set the same symmetric y-limits on every panel of the figure: fixed at
% +-pearsonYLim for standardised residuals, otherwise the largest bar
if resType == "pearson"
    yMax = pearsonYLim;
else
    yLims = cell2mat(get(resAxes(:), 'YLim'));
    yMax  = max(abs(yLims(:)));
end
set(resAxes(:), 'YLim', [-yMax yMax]);

% Figure title, left off for the manuscript version (the caption carries it)
% if resType == "pearson"
%     title(tl, "Residuals between weighted data and fitted distributions, standardised at n_{eff} " + ...
%               "(hollow bars: expected count < " + minExpCount + "; dotted: " + char(177) + seLines + " SE)", 'FontSize', 9)
% else
%     title(tl, "Residuals between weighted data and fitted distributions (observed - expected fraction per bin)", 'FontSize', 9)
% end

if saveFig
    exportgraphics(gcf, figName, 'Resolution', 300)
    disp("Saved " + figName)
end

%% KS and AD distances for every method and family
numRows    = numel(methods).*numel(families);

% Initialise the table columns
methodCol  = strings(numRows, 1);
familyCol  = strings(numRows, 1);
numEffCol  = NaN(numRows, 1);
KS_D       = NaN(numRows, 1);
KS_D_lower = NaN(numRows, 1);
KS_D_body  = NaN(numRows, 1);
KS_D_upper = NaN(numRows, 1);
AD_A2      = NaN(numRows, 1);
BICeffCol  = NaN(numRows, 1);

iRow = 0;
for i = 1:numel(methods)
    % Get the weighted data, sorted, and its effective sample size
    methodFit = d.S1.(methods(i));
    weightedC = methodFit.weightedC;
    if iscell(weightedC), weightedC = cell2mat(weightedC(:)); end
    weightedC = weightedC(isfinite(weightedC) & weightedC > 0);
    sortedC   = sort(weightedC(:));
    numRep    = numel(sortedC);
    numEff    = round(mean(methodFit.numCpairs));

    % Set up the weighted empirical CDF just before and just after each value
    ecdfBefore = ((1:numRep)' - 1)./numRep;
    ecdfAfter  = (1:numRep)'./numRep;

    % Find which sorted values fall in each quantile segment
    segCuts   = quantile(sortedC, segPct);
    lowerLog  = sortedC <  segCuts(1);
    bodyLog   = sortedC >= segCuts(1) & sortedC <= segCuts(2);
    upperLog  = sortedC >  segCuts(2);

    for j = 1:numel(families)
        iRow = iRow + 1;
        fitInfo = methodFit.(families(j)).fitInfo;

        % Calculate the KS distance: largest gap between empirical and fitted CDF,
        % checked on both sides of every jump in the empirical CDF
        cdfAtData = fittedCDF(fitInfo, families(j), sortedC);
        gapAbove  = ecdfAfter  - cdfAtData;
        gapBelow  = cdfAtData  - ecdfBefore;
        gaps      = max(gapAbove, gapBelow);
        KS_D(iRow)       = max(gaps);
        KS_D_lower(iRow) = max(gaps(lowerLog));
        KS_D_body(iRow)  = max(gaps(bodyLog));
        KS_D_upper(iRow) = max(gaps(upperLog));

        % Calculate the AD statistic by summing areas over a fine ln(nSR) grid:
        % A2 = numEff * integral of (Femp - F)^2 / (F (1 - F)) dF
        lnGrid    = -8:0.005:8;
        cdfGrid   = fittedCDF(fitInfo, families(j), exp(lnGrid));
        cdfUpper  = cdfGrid(2:end);                                     % Fitted CDF at the upper edge of each grid step
        ecdfUpper = cumsum(histcounts(log(sortedC), [-Inf, lnGrid(2:end)]))./numRep; % Empirical CDF at the same edges
        dF        = diff(cdfGrid);                                      % Fitted probability mass in each grid step
        usableLog = cdfUpper > 1e-10 & cdfUpper < 1 - 1e-10;            % Avoid dividing by zero far out in the tails
        integrand = (ecdfUpper - cdfUpper).^2 ./ (cdfUpper.*(1 - cdfUpper));
        AD_A2(iRow) = numEff.*sum(integrand(usableLog).*dF(usableLog));

        % Store the labels and BICeff alongside
        methodCol(iRow) = methodLbl(i);
        familyCol(iRow) = families(j);
        numEffCol(iRow) = numEff;
        BICeffCol(iRow) = fitInfo.BICtaeheefix;
    end
end

distanceT = table(methodCol, familyCol, numEffCol, KS_D, KS_D_lower, KS_D_body, KS_D_upper, AD_A2, BICeffCol, ...
    'VariableNames', ["method", "family", "numEff", "KS_D", "KS_D_below_p5", "KS_D_body", "KS_D_above_p95", "AD_A2", "BICeff"]);
disp("KS distance D (overall and within each quantile segment) and AD statistic A2 at numEff.")
disp("Descriptive distances only - no p-values, see header. Smaller is closer.")
disp(distanceT)

if saveTable
    writetable(distanceT, tableFile)
    disp("Saved " + tableFile)
end

%% Make the paper tables: methods as rows, families as columns
numPaperRows = numel(paperRows);
numFams      = numel(families);

% Initialise the values of the two tables
KS_paper = NaN(numPaperRows, numFams);
AD_paper = NaN(numPaperRows, numFams);

% Fill in each cell from the long table above
for i = 1:numPaperRows
    for j = 1:numFams
        rowLog = distanceT.method == paperRows(i) & distanceT.family == families(j);
        KS_paper(i,j) = distanceT.KS_D(rowLog);
        AD_paper(i,j) = distanceT.AD_A2(rowLog);
    end
end

% Label the rows with the method and the columns with the family
KS_paperT = array2table(KS_paper, 'RowNames', cellstr(paperRows), ...
    'VariableNames', cellstr(paperFamLbl));
AD_paperT = array2table(AD_paper, 'RowNames', cellstr(paperRows), ...
    'VariableNames', cellstr(paperFamLbl));
disp("Paper table: KS distance D")
disp(KS_paperT)
disp("Paper table: AD statistic A2 at numEff")
disp(AD_paperT)

if saveTable
    writetable(KS_paperT, KStableFile, 'WriteRowNames', true)
    writetable(AD_paperT, ADtableFile, 'WriteRowNames', true)
    disp("Saved " + KStableFile + " and " + ADtableFile)
end
