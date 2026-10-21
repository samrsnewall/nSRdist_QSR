function[modenSRinfo, mediannSRinfo, nSRcounts, LinNSRcounts, meanSR] = nSRBchron(corename,LabIDs, incDepths, excLabIDs, excDepths, dataLoc, S)
% This function calculates nSR for one core from a Bchron age-depth model,
% first running Bchron through Rscript if this core has no outputs in
% S.BchronFolderName (or S.BchronReDo is true). To run Bchron it reads the
% core's radiocarbon dates, filters them as for the RSR methods (if
% S.BchronFilter) and writes the Bchron input file, unless that file already
% exists and S.BchronReDo is false. It then calculates the nSR of every pair
% of dated depths from the mode (BMode) and median (BMedian) Bchron age at
% each depth, leaving out dates Bchron rejected in more than
% S.BchronReversalCriteria of its runs; from each Bchron run in turn (BSamp),
% leaving out the dates rejected in that run; and, if S.Lin2014Method, with
% the Lin2014 code (LinNSR). Raw 14C ages and errors are in kyr (Bchron is
% given them in yr); Bchron ages are in yr; depths in cm; SR in cm/kyr (the
% sed rates inside the Lin2014 code are cm/yr, its ratios are unit-free).
%
% modenSRinfo, mediannSRinfo, nSRcounts and LinNSRcounts are 3-row nSR
% matrices (nSR, dx_cm, dt_yr), each run starting with a column holding NaN,
% the run's first depth and its first age; nSRcounts holds every Bchron run
% side by side. meanSR is the mean SR (cm/kyr) from the median ages. If
% filtering rejects the core, the nSR outputs are empty and meanSR is NaN.
%
% S fields used: sandboxPath, BchronFolderName, BchronReDo, BchronFilter,
% BchronOutlier, BchronCalCurve, RscriptPath, DeltaRError, BchronDepthSpacing,
% BchronReversalCriteria, minNumberOfAges, normWithRunAve, Lin2014Method

%% Set up the Bchron folders
% Get the name of the Bchron folder to use
BchronFolder = S.BchronFolderName;

% Set up the path to the folder the Bchron data are in
BchronFolderPath = fullfile(S.sandboxPath, "BchronFolders", BchronFolder);

% If the Bchron folder doesn't yet exist, make it
if ~isfolder(BchronFolderPath)
    mkdir(BchronFolderPath)
end

% Set up the path to this core's Bchron outputs
coreDir = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Outputs",corename);

%% Run Bchron for this core if needed
% If the Bchron output folder for this core doesn't exist (or S.BchronReDo is
% true), run Bchronology
if ~isfolder(coreDir) || S.BchronReDo

    % If the Bchron "Outputs" folder doesn't exist, create it
    outputsFolder = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Outputs");
    if ~isfolder(outputsFolder)
        mkdir(outputsFolder)
    end

    % If the Bchron "Inputs" folder doesn't exist, create it
    inputsFolder = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Inputs");
    if ~isfolder(inputsFolder)
        mkdir(inputsFolder)
    end

    % Set up the name of the Bchron input file
    c14input = fullfile(S.sandboxPath, "BchronFolders", BchronFolder, "Inputs", ...
        corename + "_radiocarbon.txt");

    % --- Write the Bchron input file (if it doesn't exist, or if S.BchronReDo)
    if ~isfile(c14input) || S.BchronReDo
        % Get the radiocarbon data from our sources (either World Atlas or
        % Lin2014 dataset)
        if dataLoc == "WA"
            [age_kyr, depth_cm, ageErr_kyr, label] = getDataWA(corename, S);
        elseif dataLoc == "Lin2014"
            [age_kyr, depth_cm, ageErr_kyr, label] = getDatatxt(corename, S);
        end

        % If desired, filter the ages as for the RSR methods
        if S.BchronFilter
            [age_kyr, depth_cm, ageErr_kyr, label, noIncludeInfo, tooFewAges, ~, ~, ~] = ...
                filtering(age_kyr, depth_cm, ageErr_kyr, label, LabIDs, incDepths, excLabIDs, ...
                excDepths, corename, S);

            % Stop here if filtering of ages leads to core rejection
            if noIncludeInfo == 1 || tooFewAges == 1
                modenSRinfo   = [];
                mediannSRinfo = [];
                nSRcounts     = [];
                LinNSRcounts  = [];
                meanSR        = NaN;
                return
            end
        end

        % Convert the ages and errors to yr, as Bchron reads them
        numAges   = length(age_kyr);
        age_yr    = age_kyr.*1000;
        ageErr_yr = ageErr_kyr.*1000;

        % Set up the other columns of the Bchron input: thickness, the outlier
        % probabilities and the date type (runBchron.R passes Outlier1 to
        % Bchron and does not use Outlier2 or Type)
        outlierLikelihood = S.BchronOutlier;
        thickness_cm      = zeros(numAges,1);
        outlier1Prob      = ones(numAges,1).*outlierLikelihood;
        outlier2Prob      = ones(numAges,1).*0.001;
        dateType          = ones(numAges,1);

        % Write a new tab delimited file in the correct format for Bchron to
        % read, using the radiocarbon data we have
        c14inputT = table(label, age_yr, ageErr_yr, depth_cm, thickness_cm, outlier1Prob, ...
            outlier2Prob, dateType, 'VariableNames', ...
            ["id", "Age", "Error", "Depth", "Thickness", "Outlier1", "Outlier2", "Type"]);
        writetable(c14inputT, c14input, "Delimiter", 'tab')

    end

    % --- Run Bchron through Rscript
    % Get the calibration curve to use
    calCurve = S.BchronCalCurve;

    % Run Bchronology by sending a line to the terminal (which allows us to
    % run Rscript within MATLAB), with the name of the folder to save it in,
    % the name of the core, the calibration curve to be used, the reservoir
    % age error, and the depth spacing to use
    disp("Running Bchronology for core " + corename)
    RscriptCommand = S.RscriptPath + " " + fullfile(S.sandboxPath, "Functions", "runBchron.R") ...
        + " " + S.sandboxPath + " " + fullfile("BchronFolders", BchronFolder) ...
        + " " + corename + " " + calCurve + " " + S.DeltaRError + " " + S.BchronDepthSpacing;
    [status, Routput] = system(RscriptCommand);

    % Check whether the command ran smoothly or not
    if status == 0
        disp("|============================================================| 100%")
    else
        disp("R output below:")
        disp(Routput)
        error("Problem in Bchron")
    end
end

%% Load the Bchron outputs
% Load the radiocarbon data Bchron was given
rData = readtable(fullfile(coreDir, "inputData.txt"));

% Load the individual run data: the Bchron ages (yr) at every predicted
% position (thetaData) and whether each date was flagged as an outlier
% (phiData, one column per date), one row per run
thetaData = readmatrix(fullfile(coreDir, "theta.csv"), "NumHeaderLines",1);
phiData   = readmatrix(fullfile(coreDir, "phi.csv"), "NumHeaderLines", 1);

% Round the depth of each radiocarbon date to the S.BchronDepthSpacing used
% for the predicted positions
depthR_r2BDS = round(rData.Depth./S.BchronDepthSpacing).*S.BchronDepthSpacing;

% Set up all the depths Bchron estimated ages at
predictPositions = depthR_r2BDS(1):S.BchronDepthSpacing:depthR_r2BDS(end);

% Find the predicted positions that are radiocarbon dated, and get the
% Bchron ages (yr) at those depths
predictPositionsRLog = ismembertol(predictPositions,depthR_r2BDS);
predictPositionsR    = predictPositions(predictPositionsRLog);
thetaDataR           = thetaData(:, predictPositionsRLog);

%% Find the mode and median Bchron age at each dated depth
% Get mode and median from thetaDataR (which holds the 1000 ages used in the
% 1000 Bchron iterations for each radiocarbon dated depth - it's different
% than the radiocarbon calibration)

%%% NOTE THAT THE MODE calculated using MATLAB is often different to the
%%% MODE using Bchron. This is possible when the dataset is multimodal
%%% (which the thetaData typically are) - MATLAB will always choose the
%%% lowest mode, whilst R will (depending on your function, because there
%%% is no base R function to get the mode) get whichever Mode occurs first
%%% in the dataset. The mode is a bad choice for finding a single choice in
%%% a numerical dataset because it depends on the resolution of the data
%%% (Bchron output is likelihood per year) and there may be multiple modes.
%%%
modeMAT   = mode(thetaDataR);
medianMAT = median(thetaDataR);

% Find the percentage of runs in which each date was rejected (dividing by
% 10 assumes Bchron returned 1000 runs)
pctRejMAT = sum(phiData)./10;

% DEBUG: Test if there are any multi-modal Bchron estimates
% [~,~,test] = mode(thetaData);
% if any(cellfun(@length, test) > 1)
%     disp("Core " + corename + " has a multimodal Bchronology estimate")
% end

% Get relevant information for calculating nSR with mode of ages
modeAge_yr   = modeMAT';
medianAge_yr = medianMAT';
pctReject    = pctRejMAT';

%% Remove the ages Bchron rejected
% Remove rejected ages if desired (Lin2014 did not do this)

% Remove the age if the age was rejected more than S.BchronReversalCriteria
% of the time. Note, even in a normal looking radiocarbon profile, good ages
% are rejected up to 5% of the time
removeRejectedAges = true;
if removeRejectedAges
    % Find which ages were rejected more than S.BchronReversalCriteria of
    % the time
    rejectLog = pctReject > S.BchronReversalCriteria*100;

    % Find out what ages were accepted (using their depth info to index
    % them), and what depths had estimated ages we want to keep (unique and
    % accepted depths)
    acceptedDepths = depthR_r2BDS(~rejectLog);
    depthsUsed     = unique(acceptedDepths);

    % Keep the ages at the accepted depths (predictPositionsR holds the
    % unique dated depths, one per estimated age, so this logical is the
    % length of the ages)
    keepAgeLog       = ismembertol(predictPositionsR, depthsUsed);
    modeAgeUsed_yr   = modeAge_yr(keepAgeLog);
    medianAgeUsed_yr = medianAge_yr(keepAgeLog);
else
    depthsUsed       = depthR_r2BDS;
    modeAgeUsed_yr   = modeAge_yr;
    medianAgeUsed_yr = medianAge_yr;
end

% Handle doubly-dated depths. First option is to just find the unique
% depths, because each unique depth has a unique and singular modeAge_yr value
[depthsUsed, uniqueIdx, ~] = unique(depthsUsed);
modeAgeUsed_yr   = modeAgeUsed_yr(uniqueIdx);
medianAgeUsed_yr = medianAgeUsed_yr(uniqueIdx);

% Highlight a problem if fewer than S.minNumberOfAges ages are used
% TODO: this only warns; oneCoreScenarios treats S.minNumberOfAges as a hard cutoff
if length(depthsUsed) < S.minNumberOfAges
    warning(corename + " has less than " + num2str(S.minNumberOfAges) + " accepted ages after Bchron")
end

%% Calculate nSR with mode ages (BMode)
% Calculate the mean SR in cm/kyr (depths in cm, Bchron ages in yr, so
% multiply by 1000)
meanSR = (depthsUsed(end)-depthsUsed(1))./(modeAgeUsed_yr(end) - modeAgeUsed_yr(1)).*1000;

% Calculate SRs (cm/kyr), nSRs, dx_cm and dt_yr
SRs   = diff(depthsUsed)./diff(modeAgeUsed_yr).*1000;
nSRs  = SRs./meanSR;
dx_cm = diff(depthsUsed);
dt_yr = diff(modeAgeUsed_yr);

% Store all nSR info in one matrix, the standard set up I use
modenSRinfo = [NaN, nSRs'; depthsUsed(1), dx_cm'; modeAgeUsed_yr(1), dt_yr'];

%% Calculate nSR with median ages (BMedian)
% Calculate the mean SR, SRs, nSRs and dt_yr as for the mode ages (dx_cm is
% the same). This meanSR is the one returned, and the one BSamp normalises
% by when S.normWithRunAve is false
meanSR = (depthsUsed(end)-depthsUsed(1))./(medianAgeUsed_yr(end) - medianAgeUsed_yr(1)).*1000;
SRs    = diff(depthsUsed)./diff(medianAgeUsed_yr).*1000;
nSRs   = SRs./meanSR;
dt_yr  = diff(medianAgeUsed_yr);

% Store all nSR info in one matrix, the standard set up I use
mediannSRinfo = [NaN, nSRs'; depthsUsed(1), dx_cm'; medianAgeUsed_yr(1), dt_yr'];

%% Calculate nSR from each Bchron run (BSamp)
% Find out how many runs have fewer than S.minNumberOfAges ages accepted
numAcceptedAges   = size(phiData, 2)-sum(phiData,2);
tooFewAgesRunLog  = numAcceptedAges < S.minNumberOfAges;
numRunsTooFewAges = sum(tooFewAgesRunLog);

% Warn if that is more than a fifth of the runs
% FIXME: the warning text has no space before "runs"
if numRunsTooFewAges > size(thetaDataR,1)/5
    warning("Core " + corename + " has " + num2str(numRunsTooFewAges) + "runs with less than " ...
        + num2str(S.minNumberOfAges) + " ages")
end

% Go through each run, calculate the nSR info, and store it
for i = 1:size(thetaDataR,1)
    % Find out which ages were not rejected in this run (this also deals
    % with doubly dated depths, by finding unique depths and using
    % thetaDataR, which only estimates one age per depth)
    keepAgesLog = ~logical(phiData(i,:));

    % See if this run had enough ages (commented out because this reduces
    % the number of Bchron runs used to below 1000)
    % if sum(keepAgesLog)<S.minNumberOfAges
    %     continue
    % end

    % Get depths of radiocarbon ages not rejected in this run (one per
    % date), then the unique depths
    run_dateDepth_cm = depthR_r2BDS(keepAgesLog)';
    run_depth_cm     = unique(run_dateDepth_cm);

    % Get the ages at those depths for this run (one column of thetaDataR
    % per depth in predictPositionsR, so this logical is the length of the
    % columns)
    keepThetaRLog = ismembertol(predictPositionsR, run_depth_cm);
    run_age_yr    = thetaDataR(i,keepThetaRLog);

    % Calculate SRs (cm/kyr), dx_cm and dt_yr
    run_SRs   = diff(run_depth_cm)./diff(run_age_yr).*1000;
    run_dx_cm = diff(run_depth_cm);
    run_dt_yr = diff(run_age_yr);

    % Calculate nSR, normalising by this run's average SR or by meanSR
    if S.normWithRunAve
        run_aveSR = (run_depth_cm(end)-run_depth_cm(1))./(run_age_yr(end)-run_age_yr(1)).*1000;
        run_nSRs  = run_SRs./run_aveSR;
    else
        run_nSRs = run_SRs./meanSR;
    end

    % Store info in one block using the 3-row nSR matrix format
    nSRinfo = [NaN, run_nSRs; run_depth_cm(1), run_dx_cm; run_age_yr(1), run_dt_yr];
    if ~exist('nSRcounts','var')
        nSRcounts = nSRinfo;
    else
        nSRcounts = [nSRcounts, nSRinfo]; %#ok<AGROW>
    end
end

%% Calculate nSR with the Lin2014 method (LinNSR)
if S.Lin2014Method
    % ------ Following code reproduces data that Lin2014 code loaded in ------

    % Compute median posterior age at each predictPosition
    medAge_yr = median(thetaData, 1);   % Median ages (in yr) of each predictPosition

    % Interpolate to find depth at each 1 kyr step of median age

    % Find the relevant 1 kyr time steps (in yr, as medAge_yr is)
    age_targets_yr = ceil(min(medAge_yr)./1000).*1000:1000:floor(max(medAge_yr)./1000).*1000;

    % Make sure medAge_yr is unique - in theory it is, but because it is
    % estimated by sampling, there can be some duplication, very certain this
    % is not a serious problem
    [medAgeUnique_yr, medAgeUniqueIdx] = unique(medAge_yr, 'stable');

    % Find depths at those 1 kyr time steps
    depth_at_age = interp1(medAgeUnique_yr, predictPositions(medAgeUniqueIdx), age_targets_yr, ...
        'linear');

    % Read in median calibrated ages (yr) of each calibrated radiocarbon age
    % (not influenced by Bchron)
    calT      = readtable(fullfile(coreDir, "calAgesMedian.csv"));
    medCalAge = calT.pct50;

    % ------ Following code is exact replica of code used by Lin2014 ---------

    % Bchron age models frequently don't include first and last 14c dates
    % due to even 1-kyr increments of output
    % Add these ages back in
    if rData.Depth(1) < min(depth_at_age)
        depth_at_age=[rData.Depth(1), depth_at_age];
        age_targets_yr =[medCalAge(1), age_targets_yr];
    end
    if rData.Depth(end) > max(depth_at_age)
        depth_at_age=[depth_at_age, rData.Depth(end)];
        age_targets_yr=[age_targets_yr, medCalAge(end)];
    end

    % Bchron model ages for depth of each 14c measurement
    Bage_14c=interp1(depth_at_age,age_targets_yr,rData.Depth);

    % Calculate sed rates (NaN indicates duplicate 14C at a particular depth)
    % (cm/yr, as Bage_14c is in yr; the sed rate ratio is unit-free)
    cores_sr=diff(rData.Depth)./diff(Bage_14c);

    % Mean sed rate for core
    cores_meansr=(max(rData.Depth)-min(rData.Depth))/(max(Bage_14c)-min(Bage_14c));

    % Divide by mean sed rate to convert to sed rate ratio
    sed_rate_ratio=cores_sr/cores_meansr;
    ind1=find(~isnan(sed_rate_ratio));

    % Duration of sed rate
    depint=diff(rData.Depth);
    dur=diff(Bage_14c);
    ind2=find(dur~=0); % SN Note: ind1 and ind2 will be identical
    LinNSRcounts=[NaN, sed_rate_ratio(ind1)'; rData.Depth(1), depint(ind2)';  Bage_14c(1), dur(ind2)'];
else
    LinNSRcounts = [];
end
