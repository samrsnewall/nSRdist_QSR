function[nSRcountsArray] = countsCell2Array(nSRcountsCell, chooseLog)
% This function puts the nSR matrices of the chosen cores side by side in
% one 3-row nSR matrix. nSRcountsCell holds one core's nSR matrix per cell
% (nSR, dx_cm, dt_yr, each run starting with a NaN column; see README "nSR
% matrix format"), and the cells whose chooseLog entry is not 1 are skipped.
% nSR is unit-free; dx in cm; dt in yr.
%
% S fields used: none

%% Combine the chosen cores' nSR matrices
% Initialise the array to be added to (with a column of ones)
nSRcountsArray = ones(3,1);

% Concatenate all nSRcounts that are in the desired subset
for i = 1:length(nSRcountsCell)
    if chooseLog(i) == 1
        nSRcountsArray = cat(2, nSRcountsArray, nSRcountsCell{i});
    end
end

% Remove the column of ones that was used to set up the array
nSRcountsArray = nSRcountsArray(:,2:end);
end
