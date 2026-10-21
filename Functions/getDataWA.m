function[age_kyr, depth_cm, ageErr_kyr, label] = getDataWA(corename, S)
% This function reads the radiocarbon dates of one core from its NetCDF file
% in the Mulitza 2022 World Atlas (WA) dataset (<corename>.age in the Age
% folder of S.WApath), in the order they are in the file. Ages and their
% errors (the file's "Age +Error") are uncalibrated 14C kyr BP, as in the
% file; depths are converted from m to cm and rounded to 2 decimal places;
% label holds the LabIDs.
%
% S fields used: WApath

%% Find the core's WA file
% Set up the path to the core's .age file
WAfilePath = fullfile(S.WApath, "Age/", corename + ".age");

%% Read in radiocarbon data from the core
% Read in the depths (m) and convert them to cm, rounding to the nearest
% 100th of a cm to avoid floating point number problems
depth_m  = ncread(WAfilePath, "Depth");
depth_cm = depth_m.*100;
depth_cm = round(depth_cm, 2);

% Read in the ages and their errors (14C kyr BP)
age_kyr    = ncread(WAfilePath, "Age dated");
ageErr_kyr = ncread(WAfilePath, "Age +Error");

% Read in the LabIDs
label = ncread(WAfilePath, "Label");
label = string(label);
end
