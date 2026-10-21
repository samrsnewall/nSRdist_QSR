function [ageProb, calAge_yr] = multiMatcalQ(age_kyr, ageErr_kyr, dateIdx, S)
% This function calibrates the radiocarbon dates picked out by dateIdx
% (indices into age_kyr and ageErr_kyr) against Marine20, calling matcal
% once per date (or my faster copy, matcalFast, if S.matcalFast is true),
% and puts the results into ageProb, one column per date on the shared
% calibrated age grid calAge_yr. DeltaR is set to 0, so the Marine20 curve's
% built-in reservoir correction is used as it is, and only its error,
% S.DeltaRError (14C yr), is propagated through the calibration. Ages and
% errors come in as 14C kyr BP and are passed to matcal in 14C yr; calAge_yr
% is in calendar yr BP (every yr with matcal, every 3 yr with matcalFast)
% and each column of ageProb gives the probability of each age on that grid
% (each column sums to 1).
%
% S fields used: DeltaRError, matcalFast

%% Set up reservoir age correction
% Set DeltaR to zero and take its 1-sigma error from S (both in 14C yr)
resAge_yr   = 0;              % DeltaR is assumed zero; only the error is used
resError_yr = S.DeltaRError;  % 1-sigma uncertainty on DeltaR

%% Calibrate first age to get length of cal age vector
% Convert the first age and its error from kyr to yr, as matcal wants yr
age_yr    = age_kyr(dateIdx(1)) * 1000;
ageErr_yr = ageErr_kyr(dateIdx(1)) * 1000;

% Calibrate it with matcal or matcalFast (calProb has the calibrated ages in
% yr BP in column 1 and their probabilities in column 2)
if ~S.matcalFast
    [~, ~, calProb, ~] = matcal(age_yr, ageErr_yr, ...
        'Marine20', 'CalBP', 'resage', resAge_yr, 'reserr', resError_yr, 'plot', 0);
else
    [~, ~, calProb, ~] = matcalFast(age_yr, ageErr_yr, ...
        'Marine20', 'CalBP', 'resage', resAge_yr, 'reserr', resError_yr, 'plot', 0);
end

% Get the calibrated age grid from the first call only (it is identical for
% every call, since the calibration curve and axis are the same)
calAge_yr = calProb(:,1);

% Initialise the matrix to hold the probabilities of all the dates, and
% store the first date's
numCalAges   = length(calAge_yr);
numAges      = length(dateIdx);
ageProb      = zeros(numCalAges, numAges);
ageProb(:,1) = calProb(:,2);

%% Calibrate the rest of the ages
for i = 2:numAges
    % Convert this age and its error from kyr to yr
    age_yr    = age_kyr(dateIdx(i)) * 1000;
    ageErr_yr = ageErr_kyr(dateIdx(i)) * 1000;

    % Calibrate it and store its probabilities
    if ~S.matcalFast
        [~, ~, calProb, ~] = matcal(age_yr, ageErr_yr, ...
            'Marine20', 'CalBP', 'resage', resAge_yr, 'reserr', resError_yr, 'plot', 0);
    else
        [~, ~, calProb, ~] = matcalFast(age_yr, ageErr_yr, ...
            'Marine20', 'CalBP', 'resage', resAge_yr, 'reserr', resError_yr, 'plot', 0);
    end
    ageProb(:,i) = calProb(:,2);
end
end
