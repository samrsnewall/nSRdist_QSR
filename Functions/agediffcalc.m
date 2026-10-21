function [dtVals_yr, dtProbs, IDpairs, agediffV] = agediffcalc(ageProb, calAge_yr, ...
    numDatePairs, label, IDpairs, agediffV, S)
% This function works out the PDF of the age difference of each pair of
% neighbouring dates in a scenario (the second date's age minus the
% first's) from their calibrated age PDFs. It takes every combination of an
% age of one date and an age of the other, with the product of their
% probabilities, and sums the probabilities of combinations with the same
% age difference. To save time, ages with a probability of S.pdfMinVal or
% less are left out (assumed negligible).
%
% ageProb holds the calibrated age PDFs of the scenario's dates (one column
% per date, on the calAge_yr grid) and label their LabIDs; ages and age
% differences are in yr. dtVals_yr and dtProbs hold, for each pair, the
% possible age differences and their probabilities. Working out a pair is
% slow, so IDpairs (a string array of pair IDs, "LabID1,LabID2") and
% agediffV (a cell holding each pair's [dtVals_yr, dtProbs]) cache the
% pairs already worked out: a pair in the cache is taken from it, and a new
% pair is added to the end of both.
%
% S fields used: pdfMinVal

%% Initialise the age-difference PDFs
% One cell per pair
dtVals_yr = cell(1, numDatePairs);
dtProbs   = cell(1, numDatePairs);

%% Leave out the ages with negligible probability
% Set each calibrated PDF to NaN where it is S.pdfMinVal or less, to cut
% down the number of combinations of ages
lowProbLog              = ageProb(:,:) <= S.pdfMinVal;
ageProbNaNs             = ageProb;
ageProbNaNs(lowProbLog) = NaN;

%% Work out (or take from the cache) the age-difference PDF of each pair
for iPair = 1:numDatePairs

    % Check whether this pair of LabIDs is in the cache
    labelPair      = label(iPair) + "," + label(iPair+1);
    pairInCacheLog = ismember(IDpairs, labelPair);

    if sum(pairInCacheLog) ~= 0
        % --- Take the pair's PDF from the cache
        cachedPDF        = agediffV{pairInCacheLog};
        dtVals_yr{iPair} = cachedPDF(:,1);
        dtProbs{iPair}   = cachedPDF(:,2);

    else
        % --- Work out the pair's PDF
        % Find the ages of each date whose probability is above S.pdfMinVal
        % (not NaN in ageProbNaNs), and their probabilities
        firstAgeIdx   = find(~isnan(ageProbNaNs(:,iPair)));
        firstAges_yr  = calAge_yr(firstAgeIdx);
        firstProbs    = ageProb(firstAgeIdx, iPair);

        secondAgeIdx  = find(~isnan(ageProbNaNs(:,iPair+1)));
        secondAges_yr = calAge_yr(secondAgeIdx);
        secondProbs   = ageProb(secondAgeIdx, iPair+1);

        % Calculate the age difference of every combination of the two
        % dates' ages, and its probability (the product of the two ages'
        % probabilities)
        dt_yr      = secondAges_yr' - firstAges_yr;
        jointProbs = firstProbs .* secondProbs';

        % Sum the probabilities of the combinations with the same age
        % difference
        [dtVals_yr{iPair}, ~, uniqueDtIdx] = unique(dt_yr);
        dtProbs{iPair} = accumarray(uniqueDtIdx, jointProbs(:));

        % Add the pair to the end of the cache
        IDpairs  = [IDpairs;  labelPair];                             %#ok<AGROW>
        agediffV = [agediffV; {[dtVals_yr{iPair}, dtProbs{iPair}]}];  %#ok<AGROW>
    end
end

end
