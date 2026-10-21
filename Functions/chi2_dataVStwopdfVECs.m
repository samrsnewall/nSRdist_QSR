function [h, p, chiStat] = chi2_dataVStwopdfVECs(data, desiredSum, binN, pdfVs, fitS)
% chi2_dataVStwopdfVECs  Chi-squared goodness-of-fit test of data against
%                        one or more candidate PDFs.
%
% Bins the observed data into binN bins, downscales the bin counts to
% desiredSum (to correct for weighted replication), and computes expected
% counts from the CDF of each candidate PDF. Optionally adjusts bin edges
% to ensure minimum expected counts (via iterateBinSizeNVec), then runs
% MATLAB's chi2gof for each PDF. A diagnostic figure can be produced
% showing the observed and expected histograms alongside the p-value and
% chi-squared statistic.
%
% INPUTS
%   data       - (numeric vector) Observed data values; may be
%                log-transformed before passing (e.g. log(nSR))
%   desiredSum - (scalar) True sample size used to downscale counts; set
%                equal to numel(data) if no replication has been applied
%   binN       - (integer) Initial number of histogram bins; the function
%                may adjust this if fitS.enforceBinSizeLimits is true
%   pdfVs      - (cell array, P×1) Candidate PDFs; each cell is a struct
%                with fields:
%                  .x         (numeric vector) PDF evaluation points
%                  .px        (numeric vector) Probability density values
%                  .numParams (scalar) Number of free parameters (used to
%                             adjust degrees of freedom in chi2gof)
%                  .pdfName   (string) Display name for plot titles
%   fitS       - (struct) Settings struct. Fields used:
%                  .enforceBinSizeLimits (logical) If true, call
%                     iterateBinSizeNVec to adjust bins until each bin has
%                     an expected count of at least 5
%                  .dispChi2 (logical) If true, produce a diagnostic figure
%
% OUTPUTS
%   h       - (1 × P logical vector) Rejection decision for each PDF:
%             1 = reject H₀ (data are not drawn from that PDF) at α = 0.05
%   p       - (1 × P numeric vector) p-value for each test
%   chiStat - (1 × P cell array) chi2gof output structs, one per PDF,
%             each containing fields: chi2stat, df, edges, O, E
%
% Note: this function calls iterateBinSizeNVec, which must be present on
% the MATLAB path.
%
% See also: calcBinCounts, ARfitdists, IRfitdists

%Define some initial bins
binEdges    = linspace(min(data), max(data), binN+1);
interiorEdges = binEdges(2:end-1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;

%%% get cdfs of each pdf
for i = 1:size(pdfVs,1)
    cumsumpdfi = cumsum(pdfVs{i}.px);
    pdfVs{i}.cdf_x = cumsumpdfi./max(cumsumpdfi);
end

%Calculate observed and expected bin counts with given binEdges
[obsCountsDown, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);

fitS.enforceBinSizeLimits = true;
if fitS.enforceBinSizeLimits
    [binEdges, expCounts, obsCountsDown, binCenters] = iterateBinSizeNVec(data, pdfVs, binN, desiredSum, fitS);
end

%Perform chi2gof
for i = 1:size(pdfVs, 1)
    [h(i), p(i), chiStat{i}] = chi2gof(binCenters, 'Frequency', obsCountsDown, 'Edges', binEdges,'Expected', expCounts(i,:), 'Emin', 5, 'NParams', pdfVs{i}.numParams);
end

%If chistat has changed binedges, use those
if numel(chiStat{1}.edges) ~= numel(binEdges)
    disp('chi2gof function modified bins')
    binEdges = chiStat{1}.edges;
    binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;
    %Find observed and expected counts
    [obsCountsDown, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);
end

%Plot figure, with function, observed histogram, expected histogram, pvalue
%and chistat value
if fitS.dispChi2 == true

    %loop for figure
    figure;
    for i = 1:size(pdfVs, 1)
        subplot(size(pdfVs, 1), 1, i)
        hold on
        %Plot histograms
        %yyaxis("left")
        %Set x limits
        binEdgesDiff = (binEdges(3) - binEdges(2));
        xlim([binEdges(2)-binEdgesDiff binEdges(end-1)+binEdgesDiff])

        %Plot observed counts
        hObs = histogram('BinCounts',obsCountsDown(2:end-1), 'BinEdges', binEdges(2:end-1), 'FaceColor', 'k', "DisplayName", "Observed Counts");
        plot([-100 binEdges(2)], [obsCountsDown(1) obsCountsDown(1)], 'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")
        plot([100 binEdges(end-1)], [obsCountsDown(end) obsCountsDown(end)], 'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")

        %Plot expected counts
        hExp = histogram('BinCounts', expCounts(i,2:end-1), 'BinEdges', binEdges(2:end-1), 'FaceColor', 'r', "DisplayName","Expected Counts");
        plot([-100 binEdges(2)], [expCounts(1) expCounts(1)], 'Color', [0.9 0 0], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")
        plot([100 binEdges(end-1)], [expCounts(end) expCounts(end)], 'Color', [0.9 0 0], 'LineStyle', '-', 'LineWidth', 2, "HandleVisibility", "off")
        if numel(data) == desiredSum
            ylabel("Counts")
        else
            ylabel("Downweighted Counts")
        end

        if fitS.enforceBinSizeLimits == false
            %Plot function
            yyaxis("right")
            lTested = plot(pdfVs{i}.x, pdfVs{i}.px, 'k', 'LineWidth', 2, "DisplayName", "Tested PDF"); %#ok<*NASGU>
            ylabel("PDF")
        end

        %Add pval and chi2stat to legend
        lchi2stat = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "chi2stat = "  + num2str(chiStat{i}.chi2stat, 3));
        lpval   = plot(nan, nan, 'LineStyle', 'none', 'DisplayName', "p = " + num2str(p(i), 3) + "; h = " + num2str(h(i), 1));
        %lhnought= plot(nan, nan, 'LineStyle', 'none', "Marker", "none", 'DisplayName', "h = " + num2str(h(i), 1));
        
        legend()
        % xlim([-3 3])
        title(pdfVs{i}.pdfName)
    end
    xlabel("log(nSR)") %Note, it's not always log nSR...
end
end