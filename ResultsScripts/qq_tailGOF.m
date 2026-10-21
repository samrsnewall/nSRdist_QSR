%% qq_tailGOF — Q-Q plots of fitted nSR distributions (Reviewer 1, point 2)
%
% Response to the request for tail-focused goodness-of-fit diagnostics.
% For each nSR estimation method (rows) and each fitted family (columns),
% plots fitted quantiles (x) against weighted empirical quantiles (y) in
% ln(nSR) space, with a 1:1 line, an approximate pointwise 95 % envelope,
% and markers at the 1st / 5th / 95th / 99th percentiles of the fit.
%
% One figure, rows = methods. The dt-restricted methods (BMedian, BSamp,
% RSR500, RSR1000) share a common, tight axis range. RSR0 has no minimum dt,
% so its upper tail is far heavier; sharing its range would push the axes
% out to about +-7 and squash the other rows, so the RSR0 row gets its own
% axis limits (ownLimMethods).
%
% Why ln(nSR): nSR spans orders of magnitude, and on a linear scale the
% upper tail dominates the panel while the lower tail (slow sedimentation)
% is squashed against zero. Log space treats the two tails symmetrically;
% for the LN family it is the ordinary normal Q-Q plot.
%
% Reading a panel: points above the 1:1 line at the right-hand end mean the
% data's upper tail is heavier than the fit (empirical extreme quantiles are
% larger than the fitted ones); points below the line at the left-hand end
% mean the data's lower tail is heavier (empirical quantiles are smaller).
% Points inside the 1:1 line at both ends mean the fit's tails are too heavy.
%
% Also prints a table of empirical vs fitted percentiles (1, 5, 50, 95, 99)
% per method x family, which is the numerical version of the same figures.
%
% Run from the ResultsScripts folder (paths are relative, as in results_Newall.m).

clear;
%close all
addpath('../Functions')

%% Settings
fitFile   = "../Results/results_QSR/dataT_QSR_fit_depthweight.mat";
families  = ["LN", "MLN", "Gam", "invGam"];
familyLbl = ["Lognormal", "Mixture Lognormal", "Gamma", "Inverse Gamma"];
tailPct   = [0.01 0.05 0.95 0.99];      % percentiles marked on each panel
tablePct  = [0.01 0.05 0.50 0.95 0.99]; % percentiles reported in the table
opts.envelope = true;
opts.envLevel = 0.95;
saveFig   = true;
dateTag   = string(datetime("today", "Format", "ddMMMyy"));

% Rows of the figure, top to bottom; ownLimMethods get their own axis limits
methods       = ["BMedianLin", "BSampLin", "RSR0", "RSR500", "RSR1000"];
methodLbl     = ["BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];
ownLimMethods = "RSR0";
figName       = "../Figures/QQ_tailGOF_" + dateTag + ".png";
size_cm       = [24 27];

%% Load fit results
load(fitFile, "d")

%% Figure
ownLimLog = ismember(methods, ownLimMethods);
fig = plotQQgrid(d, methods, methodLbl, ownLimLog, families, familyLbl, tailPct, opts, size_cm);
if saveFig
    exportgraphics(fig, figName, "Resolution", 300)
    disp("Saved " + figName)
end

%% Table: empirical vs fitted percentiles (nSR, linear space), all methods
rows = cell(0, 2 + numel(tablePct));
for m = 1:numel(methods)
    wc = d.S1.(methods(m)).weightedC;
    if iscell(wc), wc = cell2mat(wc(:)); end
    wc = wc(isfinite(wc) & wc > 0);
    emp = quantile(wc, tablePct);
    rows(end+1, :) = [{methodLbl(m), "empirical"}, num2cell(emp)]; %#ok<SAGROW>
    for f = 1:numel(families)
        fit = fittedQuantile(d.S1.(methods(m)).(families(f)).fitInfo, families(f), tablePct);
        rows(end+1, :) = [{methodLbl(m), families(f)}, num2cell(fit(:)')]; %#ok<SAGROW>
    end
end
pctNames = "p" + string(round(tablePct*100));
QTable = cell2table(rows, "VariableNames", ["method", "source", pctNames]);
disp(QTable)


%% ------------------------------------------------------------------------
function fig = plotQQgrid(d, methods, methodLbl, ownLimLog, families, familyLbl, tailPct, opts, size_cm)
% One Q-Q grid: rows = methods, columns = families, axes in ln(nSR).
% Axis limits come from the quantiles actually plotted. Rows flagged in
% ownLimLog get their own limits (x and y together, so the 1:1 line stays
% diagonal); all other rows share one common range.

nM = numel(methods); nF = numel(families);

QQ = cell(nM, nF);
for m = 1:nM
    for f = 1:nF
        QQ{m, f} = qqPairs(d.S1.(methods(m)), families(f), opts);
    end
end

% Axis limits per row, rounded out to the nearest 0.5
allq   = cellfun(@(q) [log(q.qEmp); log(q.qFit)], QQ, "UniformOutput", false);
rowq   = arrayfun(@(m) vertcat(allq{m, :}), 1:nM, "UniformOutput", false);
toLim  = @(v) [floor(min(v)*2)/2, ceil(max(v)*2)/2];
rowLim = repmat(toLim(vertcat(rowq{~ownLimLog})), nM, 1);
for m = find(ownLimLog)
    rowLim(m, :) = toLim(rowq{m});
end

fig = figure("Units", "centimeters", "Position", [2 2 size_cm]);
tl  = tiledlayout(nM, nF, "TileSpacing", "compact", "Padding", "compact");

for m = 1:nM
    for f = 1:nF
        q   = QQ{m, f};
        lim = rowLim(m, :);
        ax  = nexttile; hold(ax, "on"); box(ax, "on")

        % Envelope (approximate: assumes n_eff independent draws)
        if isfield(q, "lo")
            xe = log(q.qFit);
            fill([xe; flipud(xe)], [log(q.lo); flipud(log(q.hi))], [0.85 0.85 0.85], ...
                 "EdgeColor", "none", "HandleVisibility", "off")
        end

        % 1:1 line and tail-percentile markers (at the fitted quantiles)
        plot(lim, lim, "k-", "LineWidth", 0.75, "HandleVisibility", "off")
        qt = log(fittedQuantile(d.S1.(methods(m)).(families(f)).fitInfo, families(f), tailPct));
        for t = 1:numel(qt)
            xline(qt(t), ":", "Color", [0.4 0.4 0.4], "HandleVisibility", "off")
        end

        % Q-Q points
        plot(log(q.qFit), log(q.qEmp), ".", "Color", [0.10 0.35 0.70], "MarkerSize", 5)

        xlim(lim); ylim(lim); axis square
        % text(0.04, 0.93, sprintf("BIC_{eff} = %.0f", q.BICeff), "Units", "normalized", "FontSize", 8)
        % text(0.04, 0.84, sprintf("n_{eff} = %d", q.nEff),     "Units", "normalized", "FontSize", 7)

        if m == 1,  title(familyLbl(f), "FontWeight", "normal"); end
        if f == 1,  ylabel(methodLbl(m) + newline + "empirical ln(nSR)" + newline + sprintf("n_{eff} = %d", q.nEff)); end
        if m == nM, xlabel("fitted ln(nSR)"); end
    end
end

% title(tl, "Q-Q plots of fitted nSR distributions (depth-weighted; dotted lines: 1, 5, 95, 99th fitted percentiles; " + ...
%           "grey: approx. 95 % pointwise envelope at n_{eff})", "FontSize", 9)
end
