%% quantileSegmentFits — quantile-separated fit curves (Reviewer 1, point 2)
%
% The reviewer asked for "quantile-separated fit curves for all four
% families" to quantify "segment-wise error". This script takes the usual
% histogram-plus-fitted-PDF figure and splits it into three quantile
% segments of the weighted empirical distribution — below the 5th
% percentile, the 5th–95th body, and above the 95th — drawing each segment
% on its own y-scale so tail misfit is visible instead of being a few
% pixels high under the body.
%
% Everything is in ln(nSR). The histogram is the depth-weighted empirical
% density (weightedC, the replicated vector that was fitted), normalised to
% integrate to 1 over ln(nSR). The four curves are the stored ln-space PDFs
% (<method>.<family>.lnSR), the same ones plotted in the manuscript.
%
% Error bars on the histogram are approximate binomial standard errors at
% the effective sample size n_eff = numCpairs (the same n as BICeff), so
% they show how much a bar could move if you had that many independent
% pairs, not the noise of the replicated vector. They assume independence,
% as BICeff does, and should be described as approximate.
%
% Segment-wise error is tabulated as the probability mass each fitted
% family places in each segment, against the empirical 5 / 90 / 5 %. A
% family with a slow tail that is too heavy puts more than 5 % below the
% empirical 5th percentile, and so on.
%
% One figure with all five methods (BMedian, BSamp, RSR0, RSR500, RSR1000),
% one row each, as for the Q-Q plots.
%
% Run from the ResultsScripts folder (paths are relative, as in results_Newall.m).

clear; close all
addpath('../Functions')

%% Settings
fitFile   = "../Results/results_QSR/dataT_QSR_fit_depthweight.mat";
families  = ["LN", "MLN", "Gam", "invGam"];
familyLbl = ["Lognormal", "Mixture Lognormal", "Gamma", "Inverse Gamma"];
segPct    = [0.05 0.95];   % empirical percentiles that separate the three segments
xlimSeg   = {[-2.5 -1], [-1 1], [1 2.5]};   % fixed x-limits (ln nSR) for the lower-tail, body and upper-tail columns, the same in every row and both figures
% Dashed lines mark each method's empirical 5th and 95th percentile cuts, in whichever panel they fall.
showErr   = false;         % approximate binomial error bars at n_eff
noErrBarMethods = "BMedianLin";   % methods drawn without error bars even when showErr is true
sharedY   = true;          % if true, panels in the same column (segment) share a y-axis range across methods
saveFig   = true;
dateTag   = string(datetime("today", "Format", "ddMMMyy"));

%Bin Edges for Histograms
bw           = 0.1; %bin width
logBinEdges  = -5:bw:5;          %lognSR
logBinCenters = (logBinEdges(1:end-1)+logBinEdges(2:end)).*0.5;

% Colours as in the manuscript figures (results_QSR.m): the colour-blind
% friendly scheme of Petroff (2021), LN red, MLN orange, Gam purple, invGam
% blue, in the order of families
famCol = [228  37  54; ...
          248 156  32; ...
          150  74 139; ...
           87 144 252]./255;

figSpec(1).methods = ["BMedianLin", "BSampLin", "RSR0", "RSR500", "RSR1000"];
figSpec(1).labels  = ["BMedian", "BSamp", "RSR0", "RSR500", "RSR1000"];
figSpec(1).figName = "../Figures/QuantileSegmentFits_" + dateTag + ".png";
figSpec(1).size_cm = [19.5 20.5];

%% Load fit results
load(fitFile, "d")

%% Figures
for k = 1:numel(figSpec)
    fig = plotSegmentGrid(d, figSpec(k).methods, figSpec(k).labels, families, familyLbl, ...
                          famCol, segPct, logBinEdges, xlimSeg, showErr, noErrBarMethods, sharedY, figSpec(k).size_cm);
    if saveFig
        exportgraphics(fig, figSpec(k).figName, "Resolution", 300)
        disp("Saved " + figSpec(k).figName)
    end
end

%% Table: probability mass per segment (per cent), empirical vs fitted
allMethods = [figSpec.methods];
allLabels  = [figSpec.labels];
rows = cell(0, 5);
for m = 1:numel(allMethods)
    S  = d.S1.(allMethods(m));
    wc = S.weightedC; if iscell(wc), wc = cell2mat(wc(:)); end
    wc = wc(isfinite(wc) & wc > 0);
    cut = quantile(wc, segPct);                       % linear nSR at the two cuts
    rows(end+1, :) = {allLabels(m), "empirical", 100*segPct(1), 100*(segPct(2)-segPct(1)), 100*(1-segPct(2))}; %#ok<SAGROW>
    for f = 1:numel(families)
        Fc = fittedCDF(S.(families(f)).fitInfo, families(f), cut);
        rows(end+1, :) = {allLabels(m), families(f), 100*Fc(1), 100*(Fc(2)-Fc(1)), 100*(1-Fc(2))}; %#ok<SAGROW>
    end
end
SegTable = cell2table(rows, "VariableNames", ["method", "source", "below_p5_pct", "body_pct", "above_p95_pct"]);
disp("Probability mass in each segment (%). Segments are set by the EMPIRICAL 5th and 95th percentiles,")
disp("so the empirical row is 5 / 90 / 5 by construction; fitted rows show where each family's mass actually falls.")
disp(SegTable)


%% ------------------------------------------------------------------------
function fig = plotSegmentGrid(d, methods, methodLbl, families, familyLbl, famCol, ...
                               segPct, logBinEdges, xlimSeg, showErr, noErrBarMethods, sharedY, size_cm)
% Rows = methods; columns = lower tail | body (double width) | upper tail.
% With sharedY, each column's y-range is the largest autoscaled range among
% its rows, so the same segment is comparable across methods.

nM = numel(methods); nF = numel(families);
fig = figure("Units", "centimeters", "Position", [2 2 size_cm]);
tl  = tiledlayout(nM, 4, "TileSpacing", "compact", "Padding", "compact");
axs = gobjects(nM, 3);

for m = 1:nM
    S  = d.S1.(methods(m));
    wc = S.weightedC; if iscell(wc), wc = cell2mat(wc(:)); end
    wc = wc(isfinite(wc) & wc > 0);
    lwc  = log(wc);
    N    = numel(lwc);
    nEff = round(mean(S.numCpairs));
    q    = quantile(lwc, segPct);                     % ln(nSR) cut points

    % Decide whether this method's histogram gets error bars
    drawErrLog = showErr && ~any(methods(m) == noErrBarMethods);

    % All three panels have fixed x-limits (the same for every method) and
    % use the standard histogram bins (logBinEdges), as in all other nSR
    % histograms.
    segXlim  = xlimSeg;
    segSpan  = [1 2 1];

    for s = 1:3
        ax = nexttile([1 segSpan(s)]); hold(ax, "on"); box(ax, "on")
        axs(m, s) = ax;
        e  = logBinEdges;
        bw = diff(e);
        c  = histcounts(lwc, e);
        p  = c ./ N;                                  % weighted fraction per bin
        dens = p ./ bw;                               % density in ln(nSR)

        % Histogram as grey bars (density, integrates to 1 over ln nSR)
        histogram("BinEdges", e, "BinCounts", dens, "FaceColor", [0.82 0.82 0.82], ...
                  "EdgeColor", [0.5 0.5 0.5], "DisplayName", "weighted data")

        % Approximate binomial standard error at n_eff
        if drawErrLog
            se = sqrt(p .* (1 - p) ./ nEff) ./ bw;
            xc = e(1:end-1) + bw/2;
            errorbar(xc, dens, se, "LineStyle", "none", "Color", [0.35 0.35 0.35], ...
                     "CapSize", 2, "HandleVisibility", "off")
        end

        % Fitted ln-space PDFs, clipped to the panel's x-range
        xl = segXlim{s};
        ymax = max([0, dens(e(1:end-1) < xl(2) & e(2:end) > xl(1))]);   % tallest visible bar
        if drawErrLog
            vis  = e(1:end-1) < xl(2) & e(2:end) > xl(1);
            ymax = max([ymax, dens(vis) + se(vis)]);
        end
        for f = 1:nF
            px = S.(families(f)).lnSR.px(:); xx = S.(families(f)).lnSR.x(:);
            in = xx >= xl(1) & xx <= xl(2);
            plot(xx(in), px(in), "-", "Color", famCol(f,:), "LineWidth", 1.2, "DisplayName", familyLbl(f))
            ymax = max([ymax; px(in)]);
        end

        % Mark the empirical 5th and 95th percentile cuts in whichever
        % panel's x-range they fall in (tail or body)
        for iq = 1:numel(q)
            if q(iq) >= xl(1) && q(iq) <= xl(2)
                xline(q(iq), "--", "Color", [0.3 0.3 0.3], "HandleVisibility", "off");
            end
        end

        xlim(xl)
        ylim([0, 1.05*ymax])
        if s == 1,  ylabel(methodLbl(m) + newline + "density in ln(nSR)"); end
        if m == nM, xlabel("ln(nSR)"); end
        if s == 2 && m == 1
            lgd = legend("Orientation", "horizontal", "FontSize", 9, "Box", "off");
            lgd.Layout.Tile = "south";   % place the legend below all panels
        end
        if s == 1
            text(0.04, 0.92, sprintf("n_{eff} = %d", nEff), "Units", "normalized", "FontSize", 9)
        end
    end
end

% Share the y-range within each column if requested
if sharedY
    for s = 1:3
        ymax = max(arrayfun(@(a) a.YLim(2), axs(:, s)));
        set(axs(:, s), "YLim", [0 ymax])
    end
end

end
