function[] = corePlotCal(corename, LabIDs, incDepths, excLabIDs, excDepths, dataLoc, S)
%% Read in Radiocarbon Data
%%Read in radiocarbon data from the core
if dataLoc == "WA"
    [age, depth_cm, error, label] = getDataWA(corename, S);
elseif dataLoc == "Lin2014"
    [age, depth_cm, error, label] = getDatatxt(corename, S);
end

%% Filter for plotting
[age, depth, error, ~, ~,~, EM, ManE, NotCCR, ageGapTooHigh] = filtering(age, depth_cm, error, label, LabIDs, incDepths, excLabIDs, excDepths, corename, S);

%% Plot all radiocarbon data
% figure
% hold on
% errorbar(depth, age, error, "vertical", 'o', "color", 'k', 'DisplayName', "Accepted")
% errorbar(EM.depth, EM.age, EM.error, "vertical", '.', 'color', 'r', 'DisplayName', 'Material Not Chosen')
% errorbar(ManE.depth, ManE.age, ManE.error, "vertical", 'o', 'color', 'r', 'DisplayName', "Manually Excluded")
% errorbar(NotCCR.depth, NotCCR.age, NotCCR.error, "vertical", 'x', 'color', 'r', 'DisplayName', "Outside CC Range")
% plot(ageGapTooHigh.depth, ageGapTooHigh.age, 'r-')
% %ylim([0 age(end)*1.1])
% set(gca, 'YTickLabel',get(gca,'YTick'))
% xlabel("Depth (cm)")
% ylabel(["Radiocarbon Age","(14C kyr BP)"])
% legend("location", "southeast")
% title(corename)

resAge = 0;
resError = 200; %Same as in multiMatcal currently
medCal = zeros(size(age))';
errLowCal  = zeros(size(medCal));
errHighCal = zeros(size(medCal));
for i = 1:length(age)
    [calAgeHolder, ~,~,medCal(i)] = matcal(age(i)*1000, error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    midCal(i) = sum(p95(1:2, i)/2);

    %Get the distance from the median down to the young end of the 95%
    %range and up to its old end (errorbar draws these either side of the
    %median, so the bar covers the 95% range once)
    errLowCal(i)  = medCal(i) - p95(2,i);
    errHighCal(i) = p95(1,i) - medCal(i);
end

%Find where the difference of the median cal ages are > 5000 yrs;
gap5kyrfinder = abs(diff(medCal)) > 5000;

%Set up a line joining the two ages either side of each gap, with a NaN
%between gaps so separate gaps are not joined to each other
gap5kyrDepths = [];
gap5kyrAges   = [];
for i = 1:length(gap5kyrfinder)
    if gap5kyrfinder(i)
        gap5kyrDepths = [gap5kyrDepths, depth(i), depth(i+1), NaN]; %#ok<AGROW>
        gap5kyrAges   = [gap5kyrAges, medCal(i), medCal(i+1), NaN]; %#ok<AGROW>
    end
end



EM.medCal     = [];
EM.errLowCal  = [];
EM.errHighCal = [];
for i = 1:length(EM.age)
    %Get the median and 95% range, as for the accepted dates
    [calAgeHolder, ~,~,EM.medCal(i)] = matcal(EM.age(i)*1000, EM.error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    EM.p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    EM.errLowCal(i)  = EM.medCal(i) - EM.p95(2,i);
    EM.errHighCal(i) = EM.p95(1,i) - EM.medCal(i);
end

ManE.medCal     = [];
ManE.errLowCal  = [];
ManE.errHighCal = [];
for i = 1:length(ManE.age)
    %Get the median and 95% range, as for the accepted dates
    [calAgeHolder, ~,~,ManE.medCal(i)] = matcal(ManE.age(i)*1000, ManE.error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    ManE.p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    ManE.errLowCal(i)  = ManE.medCal(i) - ManE.p95(2,i);
    ManE.errHighCal(i) = ManE.p95(1,i) - ManE.medCal(i);
end

NotCCR.medCal     = [];
NotCCR.errLowCal  = [];
NotCCR.errHighCal = [];
for i = 1:length(NotCCR.age)
    %Get the median and 95% range, as for the accepted dates
    [calAgeHolder, ~,~,NotCCR.medCal(i)] = matcal(NotCCR.age(i)*1000, NotCCR.error(i)*1000, "Marine20", "CalBP", 'resage', resAge, 'reserr', resError, 'plot', 0);
    NotCCR.p95(:,i) = [max(calAgeHolder, [], "all"), min(calAgeHolder(1:end, 1:2), [], "all"), 999];
    NotCCR.errLowCal(i)  = NotCCR.medCal(i) - NotCCR.p95(2,i);
    NotCCR.errHighCal(i) = NotCCR.p95(1,i) - NotCCR.medCal(i);
end

 %if ~isempty(gap5kyrDepths)
    %There is an age pair that has a large gap here

    figure;
    hold on
    errorbar(depth, medCal, errLowCal, errHighCal, "Marker","none" , "color", 'k', 'DisplayName', "95% CI", "LineStyle","none")
    plot(depth, medCal, "Marker", "o", "Color","k", "DisplayName", "Median Cal Age", "LineStyle","none")
    errorbar(EM.depth, EM.medCal, EM.errLowCal, EM.errHighCal, '.', "color", 'r', 'DisplayName', "Material Not Chosen")
    errorbar(ManE.depth, ManE.medCal, ManE.errLowCal, ManE.errHighCal, 'o', "color", 'r', 'DisplayName', "Manually Excluded")
    errorbar(NotCCR.depth, NotCCR.medCal, NotCCR.errLowCal, NotCCR.errHighCal, 'x', "color", 'r', 'DisplayName', "Outside CC Range")
    %upperlim = max([medCal(end), EM.medCal(end), ManE.medCal(end), NotCCR.medCal(end)])
    if ~isempty(gap5kyrDepths)
        plot(gap5kyrDepths, gap5kyrAges, 'r-', 'DisplayName', "Median Ages > 5 kyr Apart")
    end
    ylim([0 50000])
    %set(gca, 'YTickLabel',get(gca,'YTick'))
    ytickformat('%,6.4g')
    xlabel("Depth (cm)")
    ylabel(["Calibrated Age","(cal yr BP)"])
    legend("location", "southeast")
    title(corename)
 %end

end