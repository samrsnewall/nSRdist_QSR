function[age_kyr, depth_cm, ageErr_kyr, label] = getDatatxt(corename, S)
% This function reads the radiocarbon dates of one core from its
% Bchron-format text file (the Lin2014 cores): the file in
% Lin2014Cores/BchronInput under S.sandboxPath whose name contains corename,
% with columns id, Age, Error and Depth. It sorts the dates by depth,
% shallowest first. Ages and their errors are uncalibrated 14C yr BP in the
% file and are converted to 14C kyr BP; depths are in cm in the file and are
% rounded to 2 decimal places; label holds the LabIDs (the id column).
%
% S fields used: sandboxPath

%% Find the core's text file
% Set path to folder with all txt files
txtPath = fullfile(S.sandboxPath, "Lin2014Cores", "BchronInput");

% Get names of all files
folderContents = dir(txtPath);
fileNames      = strings(length(folderContents),1);
for i = 1:length(folderContents)
    fileNames(i) = string(folderContents(i).name);
end

% Find file that has core name
% FIXME: contains matches substrings, so a core name that is part of another
% core's file name would pick up both files
coreFileLog = contains(fileNames, corename);
fileName    = fileNames(coreFileLog);

%% Read in radiocarbon data from the core
% Open file and read in data, ensuring the id column is read in as a string
importOptions = detectImportOptions(fullfile(txtPath, fileName));
importOptions = setvartype(importOptions, "id", "string");
c14T          = readtable(fullfile(txtPath, fileName), importOptions);

% Ensure data is ordered by depth, shallowest first
c14T = sortrows(c14T, "Depth");

% Pick out data of interest
age_yr    = c14T.Age;
depth_cm  = c14T.Depth;
ageErr_yr = c14T.Error;
label     = string(c14T.id);

% Convert the ages and their errors from 14C yr to 14C kyr
age_kyr    = age_yr./1000;
ageErr_kyr = ageErr_yr./1000;

% Round depths to the nearest 100th of a cm to avoid floating point number
% problems
depth_cm = round(depth_cm, 2);
end
