function[binEdges, expCounts, obsCountsDown, binCenters] = iterateBinSizeNVec(data, pdfVs, binN, desiredSum, fitS)
%Approach #2 - Start with N equally spaced bins
%Define some bins
startN = binN+1;
binEdges    = linspace(min(data), max(data), startN);

%Get expected counts for those bins
[~, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);

%Set up bin min and max
minCountNumber = fitS.chi2MinCountNum;
maxCountNumber = desiredSum/4;

%Check whether the bin counts are too high, and increase number of bins
newN = startN;
while any(expCounts(:,:) > maxCountNumber, "all")
    newN = newN+1;
    binEdges    = linspace(min(data), max(data), newN);
    [~, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);
end

%Start to decrease number of bins, until there are too few
while sum(expCounts(:,:) > maxCountNumber, "all") == 0
    %Try with less bins
    newN = newN-1;
    binEdges1    = linspace(min(data), max(data), newN);
    [~, expCounts] = calcBinCounts(data, binEdges1, desiredSum, pdfVs);
end

%Add one more bin number (highest number of bins that keeps all bins
%below maxCountNumber;
if sum(expCounts > maxCountNumber, "all") > 0
    newN = newN+1;
    binEdges    = linspace(min(data), max(data), newN);
    [~, expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);
end

%Need to check which bins are too small now
%If first bin is too small
if any(expCounts(:,1) < minCountNumber)
    %Find what a first bin size that would have at least 5 counts
    %in both pdfs
    intEdgeMin = min(binEdges(2:end-1));
    if intEdgeMin >0
        intEdgeMinEst = intEdgeMin/10;
    else
        intEdgeMinEst = intEdgeMin*10;
    end
    minVals2Try = linspace(intEdgeMinEst, max(binEdges(2:end-1)), 200);
    for i = 1:size(pdfVs,1)
        pdfi = pdfVs{i};
        expCountsMini = desiredSum.*(interp1(pdfi.x, pdfi.cdf_x, minVals2Try));
        expCountsMins(i, :) = expCountsMini;
    end
    %If any tested values are way too small for all pdfs (i.e. are NaN or 0), get rid and
    %create new tested values so the resolution of tested values is
    %better
    if all(isnan(expCountsMins(:,1)) | ismembertol(expCountsMins(:,1), 0, 1e-5))
        intEdgeMinEst = minVals2Try(find(all((~isnan(expCountsMins)&~ismembertol(expCountsMins, 0, 1e-7)), 1), 1, "first"));
        minVals2Try = linspace(intEdgeMinEst,max(binEdges(2:end-1)),  200);
        for i = 1:size(pdfVs,1)
            pdfi = pdfVs{i};
            expCountsMini = desiredSum.*(interp1(pdfi.x, pdfi.cdf_x, minVals2Try));
            expCountsMins(i, :) = expCountsMini;
        end
    end
    expCountsMinsCompat = expCountsMins > minCountNumber & expCountsMins < maxCountNumber;
    minEdgeInd = find(all(expCountsMinsCompat,1), 1);
    if sum(minEdgeInd) <1
        figure()
        yyaxis left
        hold on
        for i = 1:size(pdfVs,1)
            pdfi = pdfVs{i};
            plot(pdfi.x, pdfi.px)
        end
        error("No bin edges found that work for both distributions")
    end
    minPossBinEdge = minVals2Try(minEdgeInd);   %Find a possible bin edge that would provide at least 5 counts to minimum bin
    minBinEdge = min(binEdges(binEdges>minPossBinEdge)); %Choose even spaced bin edge that is next biggest than minimum possible bin edge
    binEdges = [binEdges(1), binEdges(binEdges >=minBinEdge)];
    [obsCountsDown,expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);
end
%If last bin is too small
if any(expCounts(:,end) < minCountNumber | isnan(expCounts(:,end)))
    %Find what last bin size would have at least minCountNumber counts
    %in all pdfs
    intEdgeMax = max(binEdges(2:end-1)); %Get the max interior bin
    % increase it to make sure it covers everything
    if intEdgeMax >0
        intEdgeMaxEst = intEdgeMax*10;
    else
        intEdgeMaxEst = intEdgeMax/10;
    end
    %Get values to test
    maxVals2Try = linspace(min(binEdges(2:end-1)), intEdgeMaxEst, 200);
    %Test them for all pdfs
    for i = 1:size(pdfVs,1)
        pdfi = pdfVs{i};
        expCountsMaxi = desiredSum.*(interp1(pdfi.x, pdfi.cdf_x, maxVals2Try));
        expCountsMaxs(i, :) = expCountsMaxi;
    end
    %If any tested values are way too large for all pdfs, get rid and
    %create new tested values so the resolution of tested values is
    %better
    if all(isnan(expCountsMaxs(:,end)) | ismembertol(expCountsMaxs(:,end), desiredSum, 1e-5))
        intEdgeMaxEst = maxVals2Try(find(all((~isnan(expCountsMaxs)&~ismembertol(expCountsMaxs, desiredSum, 1e-5)), 1), 1, "last"));
        maxVals2Try = linspace(min(binEdges(2:end-1)), intEdgeMaxEst, 200);
        for i = 1:size(pdfVs,1)
            pdfi = pdfVs{i};
            expCountsMaxi = desiredSum.*(interp1(pdfi.x, pdfi.cdf_x, maxVals2Try));
            expCountsMaxs(i, :) = expCountsMaxi;
        end
    end
    %Find which tested values fit requirements
    expCountsMaxsCompat = expCountsMaxs > desiredSum - maxCountNumber & expCountsMaxs < desiredSum - minCountNumber;
    maxEdgeInd = find(all(expCountsMaxsCompat,1), 1, 'last');
    %If there is no index, no values fit all requirements
    if sum(maxEdgeInd) <1
        figure()
        yyaxis left
        hold on
        for i = 1:size(pdfVs,1)
            pdfi = pdfVs{i};
            plot(pdfi.x, pdfi.px)
        end
        error("No bin edges found that work for both distributions")
    end
    maxPossBinEdge = maxVals2Try(maxEdgeInd);   %Find a possible bin edge that would provide at least 5 counts to maximum bin
    maxBinEdge = max(binEdges(binEdges<maxPossBinEdge)); %Choose even spaced bin edge that is next smallest than maximum possible bin edge
    binEdges = [binEdges(binEdges <= maxBinEdge), binEdges(end)];
    [obsCountsDown,expCounts] = calcBinCounts(data, binEdges, desiredSum, pdfVs);
end
interiorEdges = binEdges(2:end-1);
binCenters  = (binEdges(1:end-1)+binEdges(2:end))./2;