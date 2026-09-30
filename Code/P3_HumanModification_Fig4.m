%% Figure 4
% Calculate Figure 4 directly from the original pixel-level GHM values.

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(fileparts(codeDir), 'Statistics', 'GHM-related');
outputDir = fullfile(fileparts(codeDir), 'Analysis_Figures');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

changeFile = fullfile(dataDir, 'Change_Pixels_with_GHM.mat');
marshFile = fullfile(dataDir, 'Marsh_Pixels_with_GHM.mat');
assert(isfile(changeFile), 'Missing input file: %s', changeFile);
assert(isfile(marshFile), 'Missing input file: %s', marshFile);

edges = linspace(0, 1, 21);
centers = (edges(1:end-1) + edges(2:end)) / 2;
nBins = numel(centers);
pixelAreaHa = 0.09;

% Calculate permanent gain and loss directly from pixel-level GHM values.
S = load(changeFile, 'TidalMarshChangeRecords');
change = S.TidalMarshChangeRecords;
ghm = double(change.GHM);
category = lower(string(change.changeCategory));
use = change.year >= 1989 & change.year <= 2023 & ...
    isfinite(ghm) & ghm >= 0 & ghm <= 1;
gainGHM = ghm(use & contains(category, 'permanentgain'));
lossGHM = ghm(use & contains(category, 'permanentloss'));
assert(~isempty(gainGHM) && ~isempty(lossGHM), ...
    'No permanent gain or loss records found.');

% Calculate annual marsh area directly from all marsh-pixel GHM values.
S = load(marshFile, 'GHM_All');
marshYears = double(S.GHM_All.Year);
marshGHM = double(S.GHM_All.GHM);
yearsShown = sort(unique(marshYears(isfinite(marshYears))));
areaByYear = zeros(nBins, numel(yearsShown));
for i = 1:numel(yearsShown)
    areaByYear(:, i) = pixelAreaHa * histcounts( ...
        marshGHM(marshYears == yearsShown(i)), edges)';
end
meanMarshArea = mean(areaByYear, 2);
meanMarshArea(meanMarshArea == 0) = NaN;

changeArea = pixelAreaHa * [histcounts(gainGHM, edges)' ...
    histcounts(lossGHM, edges)'];
changePercent = 100 * changeArea ./ meanMarshArea;

% Bootstrap mapped gain and loss pixels separately.
nBoot = 400;
rng(42, 'twister');
bootArea = zeros(nBoot, nBins, 2);
for b = 1:nBoot
    sampledGain = gainGHM(randi(numel(gainGHM), numel(gainGHM), 1));
    sampledLoss = lossGHM(randi(numel(lossGHM), numel(lossGHM), 1));
    sampledArea = pixelAreaHa * [histcounts(sampledGain, edges)' ...
        histcounts(sampledLoss, edges)'];
    bootArea(b, :, :) = reshape(sampledArea, 1, nBins, 2);
end
areaLow = squeeze(prctile(bootArea, 2.5, 1));
areaHigh = squeeze(prctile(bootArea, 97.5, 1));
percentLow = 100 * areaLow ./ meanMarshArea;
percentHigh = 100 * areaHigh ./ meanMarshArea;
clear S change category ghm marshYears marshGHM use;

fontName = 'Arial';
fontSize = 18;
gainColor = [0 114 178] / 255;
lossColor = [230 159 0] / 255;
colors = [gainColor; lossColor];

fig = figure('Color', 'w', 'Position', [60 80 1450 570]);

% a. Tidal-marsh extent across GHM levels.
axA = axes(fig, 'Position', [0.063 0.16 0.262 0.68]);
hold(axA, 'on');
cmap = parula(numel(yearsShown));
for i = 1:numel(yearsShown)
    y = areaByYear(:, i);
    y(y <= 0) = NaN;
    plot(axA, centers, y, 'Color', cmap(i, :), 'LineWidth', 1.5);
end
set(axA, 'YScale', 'log');
colormap(axA, cmap);
clim(axA, [yearsShown(1)-2.5 yearsShown(end)+2.5]);
cb = colorbar(axA, 'north');
cb.AxisLocation = 'in';
cb.Ticks = yearsShown;
cb.TickLabels = string(yearsShown);
cb.Label.String = 'Year';
cb.Position = [0.0744 0.80 0.239 0.024];
axA.Position = [0.063 0.16 0.262 0.68];
ylabel(axA, 'Tidal marsh extent (ha)');

% b. Gain and loss as a percentage of marsh area in each GHM bin.
axB = axes(fig, 'Position', [0.391 0.16 0.262 0.68]);
hold(axB, 'on');
barsB = bar(axB, centers, changePercent, 'grouped', 'BarWidth', 0.9);
for i = 1:2
    barsB(i).FaceColor = colors(i, :);
    barsB(i).EdgeColor = 'none';
end
drawnow;
for i = 1:2
    errorbar(axB, barsB(i).XEndPoints, changePercent(:, i), ...
        max(0, changePercent(:, i) - percentLow(:, i)), ...
        max(0, percentHigh(:, i) - changePercent(:, i)), ...
        'k', 'LineStyle', 'none', 'CapSize', 4, ...
        'HandleVisibility', 'off');
end
set(axB, 'YLim', [0 22], 'YTick', 0:5:20);
ylabel(axB, 'Marsh change (%)');
legend(axB, barsB, {'Tidal marsh gain (%)', 'Tidal marsh loss (%)'}, ...
    'Location', 'northeast', 'Box', 'off');

% c. Absolute gain and loss area in each GHM bin.
axC = axes(fig, 'Position', [0.719 0.16 0.262 0.68]);
hold(axC, 'on');
positive = [changeArea(changeArea > 0); areaLow(areaLow > 0)];
baseValue = 10^floor(log10(min(positive)));
areaForPlot = changeArea;
areaForPlot(areaForPlot <= 0) = NaN;
barsC = bar(axC, centers, areaForPlot, 'grouped', 'BarWidth', 0.9, ...
    'BaseValue', baseValue);
for i = 1:2
    barsC(i).FaceColor = colors(i, :);
    barsC(i).EdgeColor = 'none';
end
set(axC, 'YScale', 'log', 'YLim', [1 1e4]);
drawnow;
for i = 1:2
    errorbar(axC, barsC(i).XEndPoints, areaForPlot(:, i), ...
        max(0, areaForPlot(:, i) - max(areaLow(:, i), baseValue)), ...
        max(0, areaHigh(:, i) - changeArea(:, i)), ...
        'k', 'LineStyle', 'none', 'CapSize', 4, ...
        'HandleVisibility', 'off');
end
ylabel(axC, 'Gain / loss area (ha)');
legend(axC, barsC, {'Tidal marsh gain (ha)', 'Tidal marsh loss (ha)'}, ...
    'Location', 'northwest', 'Box', 'off');

axesList = [axA axB axC];
for i = 1:3
    ax = axesList(i);
    set(ax, 'FontName', fontName, 'FontSize', fontSize, ...
        'LineWidth', 1, 'XDir', 'reverse', 'XLim', [0 1], ...
        'XTick', 0:0.2:1, 'Layer', 'top', 'YGrid', 'on', ...
        'YMinorTick', 'off', 'YMinorGrid', 'off', ...
        'XGrid', 'off', 'Box', 'on');
    xlabel(ax, 'GHM index');
    annotation(fig, 'textbox', ...
        [ax.Position(1)-0.05 0.81 0.03 0.055], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none', ...
        'Margin', 0, 'FontName', fontName, ...
        'FontSize', fontSize + 3, 'FontWeight', 'bold');
end

exportgraphics(fig, fullfile(outputDir, 'Fig4_HumanModification.png'), ...
    'Resolution', 600);
exportgraphics(fig, fullfile(outputDir, 'Fig4_HumanModification.pdf'), ...
    'ContentType', 'vector');