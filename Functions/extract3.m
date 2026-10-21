function[LabIDs, incDepths, excLabIDs, excDepths, dataLoc] = extract3(rawdataManual, chosenCoresLog, S)
% This function gets, for each core chosen by chosenCoresLog (a logical
% over the rows of the metadata table rawdataManual), the information that
% filtering needs to choose and exclude the core's radiocarbon dates, and
% where the core's dates are read from. It is the third version of this
% function (extract1 and extract2 are the earlier ones).
%
% LabIDs is "all" and incDepths is "" for every core. excLabIDs and
% excDepths are comma-separated lists (char) of the LabIDs and depths to
% exclude, made from the metadata columns chosen by the core's WAuse and
% Lin2014Keep flags and by S.removeLargeGaps and S.modifyLin2014Data (see
% each branch below); a core with nothing to exclude gets '' or []. dataLoc
% is "WA" or "Lin2014" ("" for a core with neither flag set). Depths in
% excDepths are in m, as text (filtering converts them to cm).
%
% S fields used: useLin, usePF, modifyLin2014Data, removeLargeGaps

%% Initialise the outputs
% Initialise cells
numCores  = sum(chosenCoresLog);
LabIDs    = cell(numCores, 1);
incDepths = cell(numCores, 1);
excLabIDs = cell(numCores, 1);
excDepths = cell(numCores, 1);

% Set up all LabIDs to be "all" (and incDepths to be empty)
for i = 1:numCores
    LabIDs{i}    = "all";
    incDepths{i} = "";
end

% Get the metadata of the chosen cores
rawdataUse = rawdataManual(chosenCoresLog, :);

% Initialise how to access the data (dataLoc)
dataLoc = strings(numCores,1);

%% Get the exclusions when only Lin2014 cores are used
% If we are only using Lin 2014 data, check whether we want to modify,
% get modifications and leave function
if S.useLin && ~S.usePF
    dataLoc(rawdataUse.Lin2014 == 1) = "Lin2014";
    if S.modifyLin2014Data
        Lin2014Removals = table2cell(rawdataUse(:, "Lin2014ManualExclude"));
        for i = 1:numCores
            % Join the core's Lin2014ManualExclude LabIDs into one
            % comma-separated list
            excLabIDsVec    = [string(Lin2014Removals(i))];
            excLabIDsVec    = excLabIDsVec(excLabIDsVec ~= "");
            excLabIDsJoined = "";
            for j = 1:length(excLabIDsVec)
                if j == 1
                    excLabIDsJoined = excLabIDsVec(j);
                else
                    excLabIDsJoined = excLabIDsJoined + ", " + excLabIDsVec(j);
                end
            end
            %excLabIDsJoined = [excLabIDsVec(1) + ", " + excLabIDsVec(2) + ", " ...
            %    + excLabIDsVec(3) + ", " + excLabIDsVec(4)];
            excLabIDs{i} = char(excLabIDsJoined);
        end
    end
    return
end

%% Get where each core's data are read from
% Add info about how to access the data (a core with Lin2014Keep set is
% read from the Lin2014 file even if WAuse is also set)
WAuseLog          = rawdataUse.WAuse == 1;
dataLoc(WAuseLog) = "WA";

Lin2014KeepLog          = rawdataUse.Lin2014Keep == 1;
dataLoc(Lin2014KeepLog) = "Lin2014";

%% Get the metadata columns listing the dates to exclude
% Set up excLabIDs to be a combination of all reasons to exclude LabIDs
% (depending on decision to remove large age gaps or not)
reversalIDs      = table2cell(rawdataUse(:, "ReversalIDs"));
AgeGapIDs        = table2cell(rawdataUse(:, "AgeGapIDs"));
NonPlanktonicIDs = table2cell(rawdataUse(:, "NonPlanktonicIDs"));
MiscRemovalIDs   = table2cell(rawdataUse(:, "MiscRemovalIDs"));
Lin2014Removals  = table2cell(rawdataUse(:, "Lin2014ManualExclude"));

% Set up excDepths to be a combination of all reasons to exclude depths
% (depths in m)
ReversalDepths = table2cell(rawdataUse(:, "ReversalDepths"));
AgeGapDepths   = table2cell(rawdataUse(:, "AgeGapDepths"));

%% Make each core's lists of LabIDs and depths to exclude
for i = 1:numCores
    if WAuseLog(i) == 1 && Lin2014KeepLog(i) ~= 1
        % --- WA cores (WAuse set, Lin2014Keep not set)
        % Combine the LabIDs to exclude (the age gap ones only if
        % S.removeLargeGaps)
        if S.removeLargeGaps
            excLabIDsVec = [string(reversalIDs{i}), string(AgeGapIDs{i}), ...
                string(NonPlanktonicIDs{i}), string(MiscRemovalIDs{i})];
        else
            excLabIDsVec = [string(reversalIDs{i}), string(NonPlanktonicIDs{i}), ...
                string(MiscRemovalIDs{i})];
        end

        % Join them into one comma-separated list
        excLabIDsVec    = excLabIDsVec(excLabIDsVec ~= "");
        excLabIDsJoined = "";
        for j = 1:length(excLabIDsVec)
            if j == 1
                excLabIDsJoined = excLabIDsVec(j);
            else
                excLabIDsJoined = excLabIDsJoined + ", " + excLabIDsVec(j);
            end
        end
        %excLabIDsJoined = [excLabIDsVec(1) + ", " + excLabIDsVec(2) + ", " ...
        %    + excLabIDsVec(3) + ", " + excLabIDsVec(4)];
        excLabIDs{i} = char(excLabIDsJoined);

        % Replace a NaN (no reversal depths) with an empty string
        if isnan(ReversalDepths{i})
            ReversalDepths{i} = "";
        end

        % Check settings for choice about removing large gaps or not, and
        % combine the depths to exclude
        if S.removeLargeGaps
            excDepthsVec = [string(ReversalDepths{i}), string(AgeGapDepths{i})];
        else
            excDepthsVec = [string(ReversalDepths{i})];
        end

        % Join them into one comma-separated list
        % FIXME: a cell that ends in a comma (ReversalDepths '3.78,' for
        % GeoB3304-5) gives '3.78,, 0.4', which filtering reads as NaN and
        % 0.4, so the reversal date at 378 cm is not excluded
        excDepthsVec    = excDepthsVec(excDepthsVec ~= "");
        excDepthsJoined = "";
        for j = 1:length(excDepthsVec)
            if j == 1
                excDepthsJoined = excDepthsVec(j);
            else
                excDepthsJoined = excDepthsJoined + ", " + excDepthsVec(j);
            end
        end
        excDepths{i} = char(excDepthsJoined);

    elseif WAuseLog(i) == 1 && Lin2014KeepLog(i) == 1
        % --- Cores with both WAuse and Lin2014Keep set (read from Lin2014)
        % Join the Lin2014ManualExclude LabIDs into one comma-separated list,
        % if S.modifyLin2014Data (no depths are excluded)
        if S.modifyLin2014Data
            excLabIDsVec    = [string(Lin2014Removals{i})];
            excLabIDsVec    = excLabIDsVec(excLabIDsVec ~= "");
            excLabIDsJoined = "";
            for j = 1:length(excLabIDsVec)
                if j == 1
                    excLabIDsJoined = excLabIDsVec(j);
                else
                    excLabIDsJoined = excLabIDsJoined + ", " + excLabIDsVec(j);
                end
                % FIXME: set inside the loop, so a core with no
                % Lin2014ManualExclude LabIDs keeps [] here rather than the
                % '' the other branches give (filtering treats both as empty)
                excLabIDs{i} = char(excLabIDsJoined);
            end
        end
    elseif WAuseLog(i) ~= 1 && Lin2014KeepLog(i) == 1
        % --- Lin2014 cores (Lin2014Keep set, WAuse not set)
        % Same as for the WA cores, but only if S.modifyLin2014Data
        if S.modifyLin2014Data
            % Combine the LabIDs to exclude (the age gap ones only if
            % S.removeLargeGaps)
            if S.removeLargeGaps
                excLabIDsVec = [string(reversalIDs{i}), string(AgeGapIDs{i}), ...
                    string(NonPlanktonicIDs{i}), string(MiscRemovalIDs{i})];
            else
                excLabIDsVec = [string(reversalIDs{i}), string(NonPlanktonicIDs{i}), ...
                    string(MiscRemovalIDs{i})];
            end

            % Join them into one comma-separated list
            excLabIDsVec    = excLabIDsVec(excLabIDsVec ~= "");
            excLabIDsJoined = "";
            for j = 1:length(excLabIDsVec)
                if j == 1
                    excLabIDsJoined = excLabIDsVec(j);
                else
                    excLabIDsJoined = excLabIDsJoined + ", " + excLabIDsVec(j);
                end
            end
            %excLabIDsJoined = [excLabIDsVec(1) + ", " + excLabIDsVec(2) + ", " ...
            %    + excLabIDsVec(3) + ", " + excLabIDsVec(4)];
            excLabIDs{i} = char(excLabIDsJoined);

            % Replace a NaN (no reversal depths) with an empty string
            if isnan(ReversalDepths{i})
                ReversalDepths{i} = "";
            end

            % Check settings for choice about removing large gaps or not,
            % and combine the depths to exclude
            if S.removeLargeGaps
                excDepthsVec = [string(ReversalDepths{i}), string(AgeGapDepths{i})];
            else
                excDepthsVec = [string(ReversalDepths{i})];
            end

            % Join them into one comma-separated list
            excDepthsVec    = excDepthsVec(excDepthsVec ~= "");
            excDepthsJoined = "";
            for j = 1:length(excDepthsVec)
                if j == 1
                    excDepthsJoined = excDepthsVec(j);
                else
                    excDepthsJoined = excDepthsJoined + ", " + excDepthsVec(j);
                end
            end
            excDepths{i} = char(excDepthsJoined);
        end
    end
end
