%% Extended Data Figure 4. Tidal marsh change and basin-mean human pressure
% Permanent gain and loss, 1989-2023. Observed values, without uncertainty.
% Basin-mean GHM values and histogram bins are calculated in this script.

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
figDir = fileparts(codeDir);
srcDir = fullfile(figDir, 'Statistics', 'GHM-related');
outDir = fullfile(fileparts(codeDir), 'Analysis_Figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end

changeFile = fullfile(srcDir, ...
    'Change_Pixels_with_GHM_by_Basin_1986_2024.csv');
basinGHMFile = fullfile(srcDir, ...
    'Basin_GHM_Statistics_1990_2020.csv');
basinMarshFile = fullfile(srcDir, ...
    'Annual_Area_by_Basin_1985_2024.csv');
assert(isfile(changeFile), 'Missing input: %s', changeFile);
assert(isfile(basinGHMFile), 'Missing input: %s', basinGHMFile);
assert(isfile(basinMarshFile), 'Missing input: %s', basinMarshFile);

change = readtable(changeFile, 'TextType', 'string');
basinGHM = readtable(basinGHMFile);
basinMarsh = readtable(basinMarshFile);

yearsShown = (1990:5:2020)';
edges = (0:0.025:0.6)';
centers = edges(1:end-1) + diff(edges) / 2;
nBins = numel(centers);

% Calculate annual marsh area in bins of basin-mean GHM.
validGHM = isfinite(basinGHM.basin_ID) & isfinite(basinGHM.ghm_year) & ...
    isfinite(basinGHM.human_pressure_mean) & ...
    basinGHM.human_pressure_mean >= 0 & basinGHM.human_pressure_mean <= 1;
basinGHM = basinGHM(validGHM, ...
    {'basin_ID', 'ghm_year', 'human_pressure_mean'});
[group, basinID, ghmYear] = findgroups( ...
    basinGHM.basin_ID, basinGHM.ghm_year);
ghmMean = splitapply(@mean, basinGHM.human_pressure_mean, group);
basinGHM = table(basinID, ghmYear, ghmMean, ...
    'VariableNames', {'basin_ID', 'ghm_year', 'basin_mean_GHM'});

basinMarsh.basin_ID = basinMarsh.HYBAS_ID;
basinMarsh.ghm_year = basinMarsh.year;
basinMarsh.marsh_area_ha = basinMarsh.marsh_km2 * 100;
basinMarsh = basinMarsh(ismember(basinMarsh.ghm_year, yearsShown), ...
    {'basin_ID', 'ghm_year', 'marsh_area_ha'});
areaSource = innerjoin(basinMarsh, basinGHM, ...
    'Keys', {'basin_ID', 'ghm_year'});
areaBin = discretize(areaSource.basin_mean_GHM, edges);
areaByYear = zeros(nBins, numel(yearsShown));
for i = 1:numel(yearsShown)
    use = areaSource.ghm_year == yearsShown(i) & isfinite(areaBin);
    areaByYear(:, i) = accumarray(areaBin(use), ...
        areaSource.marsh_area_ha(use), [nBins 1], @sum, 0);
end
representativeArea = mean(areaByYear, 2);

% Assign each change pixel the mean GHM of its basin in the matched GHM year.
missingPixelArea = ~isfinite(change.pixel_area_ha);
change.pixel_area_ha(missingPixelArea) = 0.09;
validChange = isfinite(change.basin_ID) & ...
    change.year >= 1989 & change.year <= 2023 & ...
    ismember(change.changeType, [-1 1]) & isfinite(change.GHMYear);
change = change(validChange, :);
changeGHM = renamevars(basinGHM, 'ghm_year', 'GHMYear');
change = innerjoin(change, changeGHM, ...
    'Keys', {'basin_ID', 'GHMYear'});
change = change(isfinite(change.basin_mean_GHM) & ...
    change.basin_mean_GHM >= 0 & change.basin_mean_GHM <= 1, :);
changeBin = discretize(change.basin_mean_GHM, edges);
gainRows = lower(change.changeCategory) == "permanentgain" & ...
    isfinite(changeBin);
lossRows = lower(change.changeCategory) == "permanentloss" & ...
    isfinite(changeBin);
gainArea = accumarray(changeBin(gainRows), ...
    change.pixel_area_ha(gainRows), [nBins 1], @sum, 0);
lossArea = accumarray(changeBin(lossRows), ...
    change.pixel_area_ha(lossRows), [nBins 1], @sum, 0);
changeArea = [gainArea lossArea];
changedPct = 100 * changeArea ./ representativeArea;
changedPct(representativeArea <= 0, :) = 0;

fontName = 'Arial';
fontSize = 18;
gainColor = [0 114 178] / 255;
lossColor = [230 159 0] / 255;
colors = [gainColor; lossColor];
fig = figure('Color', 'w', 'Position', [60 80 1450 570]);

%% a. Marsh area binned by basin-mean GHM
axA = axes(fig, 'Position', [0.063 0.16 0.262 0.68]);
hold(axA, 'on');
cmap = parula(numel(yearsShown));
for i = 1:numel(yearsShown)
    plot(axA, centers, areaByYear(:, i) / 1000, ...
        'Color', cmap(i, :), 'LineWidth', 2.5);
end
colormap(axA, cmap);
clim(axA, [yearsShown(1)-2.5 yearsShown(end)+2.5]);
cb = colorbar(axA, 'north');
cb.AxisLocation = 'in';
cb.Ticks = yearsShown;
cb.TickLabels = string(yearsShown);
cb.Label.String = 'Year';
cb.Position = [0.0744 0.80 0.239 0.024];
axA.Position = [0.063 0.16 0.262 0.68];
ylabel(axA, 'Tidal marsh extent (K ha)');
areaTop = ceil(max(areaByYear, [], 'all') / 1000 * 1.3 / 5) * 5;
set(axA, 'YScale', 'linear', 'YLim', [0 areaTop], ...
    'YTick', 0:5:areaTop-5);
axA.YRuler.Exponent = 0;

%% b. Gain and loss relative to marsh area in each basin-mean GHM bin
axB = axes(fig, 'Position', [0.391 0.16 0.262 0.68]);
hold(axB, 'on');
barsB = bar(axB, centers, changedPct, 'grouped', 'BarWidth', 0.9);
for i = 1:2
    barsB(i).FaceColor = colors(i, :);
    barsB(i).EdgeColor = 'none';
end
ylabel(axB, 'Marsh change (%)');
set(axB, 'YLim', [0 55], 'YTick', 0:10:50);
legend(axB, barsB, {'Tidal marsh gain (%)', 'Tidal marsh loss (%)'}, ...
    'Location', 'northeast', 'Box', 'off');

%% c. Gain and loss area in each basin-mean GHM bin
axC = axes(fig, 'Position', [0.719 0.16 0.262 0.68]);
hold(axC, 'on');
area = changeArea;
positiveArea = area(area > 0);
lowPower = floor(log10(min(positiveArea)));
barBase = 10^lowPower;
areaForPlot = area;
areaForPlot(areaForPlot <= 0) = NaN;
barsC = bar(axC, centers, areaForPlot, 'grouped', 'BarWidth', 0.9, ...
    'BaseValue', barBase);
for i = 1:2
    barsC(i).FaceColor = colors(i, :);
    barsC(i).EdgeColor = 'none';
end
set(axC, 'YScale', 'log', 'YLim', [barBase 1e4], ...
    'YTick', 10.^(lowPower:4));
ylabel(axC, 'Gain / loss area (ha)');
legend(axC, barsC, {'Tidal marsh gain (ha)', 'Tidal marsh loss (ha)'}, ...
    'Location', 'northwest', 'Box', 'off');

%% Labels and export
axesList = [axA axB axC];
for i = 1:3
    ax = axesList(i);
    set(ax, 'FontName', fontName, 'FontSize', fontSize, ...
        'LineWidth', 1, 'XDir', 'reverse', 'XLim', [0 0.6], ...
        'XTick', 0:0.1:0.6, 'Layer', 'top', 'YGrid', 'on', ...
        'YMinorTick', 'off', 'YMinorGrid', 'off', ...
        'XGrid', 'off', 'Box', 'on');
    xlabel(ax, 'Basin-mean GHM index');
    annotation(fig, 'textbox', [ax.Position(1)-0.05 0.81 0.03 0.055], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none', 'Margin', 0, ...
        'FontName', fontName, 'FontSize', fontSize + 3, 'FontWeight', 'bold');
end
drawnow;
exportgraphics(fig, fullfile(outDir, 'ExtendedDataFig4_BasinHumanPressure.png'), ...
    'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'ExtendedDataFig4_BasinHumanPressure.pdf'), 'ContentType', 'vector');
fprintf('Extended Data Fig. 4: %.2f ha gain, %.2f ha loss; no uncertainty intervals.\n', ...
    sum(area(:, 1)), sum(area(:, 2)));
