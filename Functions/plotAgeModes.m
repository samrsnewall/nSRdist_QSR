function[] = plotAgeModes(chosenCoresLog, interestCoresLog, ageModesCell, cores)
%This function plots the modes of the calibrated radiocarbon ages of each core on a horizontal line (where
%x value is depth of the date). This shows all the cores and all the ages
%that are being used, hence is a good simple representation of the data

% Get age modes (which were calculated in calcData.m)
ageModes = cell(size(ageModesCell));
for i = 1:length(ageModesCell)
    holder1 = ageModesCell{i};
    holder2 = NaN(0,0);
    for j = 1:length(holder1)
        holder2 = [holder2; holder1{j}];
    end
    ageModes{i} = unique(holder2);
end

%Find which indexes are of interest
int1 = chosenCoresLog + interestCoresLog;
%int2 = int1(int1 ~=0);
intIndexes = find(int1 == 2);

figure;
yTL = string(); %This will hold the core names in a format to be used for yticklabels
indexes = intIndexes;
for i = 1:length(indexes)
    %For each core, plot the mode of each calibrated radiocarbon date on a
    %horizontal line, with a 1 pt vertical offset to compare cores
    if max(diff(ageModes{indexes(i)})) < 4000
        plot(ageModes{indexes(i)}./1000, i*ones(size(ageModes{indexes(i)})), 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'Color', 'k', 'LineWidth', 1)
    else
        plot(ageModes{indexes(i)}./1000, i*ones(size(ageModes{indexes(i)})), 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'Color', [0.01, 0.01, 0.01], 'LineWidth', 1)
    end
    hold on
    xlabel("Age (kyr)")
    yTL = [yTL, cores{indexes(i)}];
end
yticks(1:length(yTL))
yticklabels(yTL(2:end))
ylim([0 length(yTL)+1])
end