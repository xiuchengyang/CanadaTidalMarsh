%% Figure 1 and Extended Data Figure 2
% Both figures are created from the same area table.

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
dataFile = fullfile(fileparts(codeDir), 'Statistics', 'Area-related', ...
    'Annual_Area_by_Coast_1988_2023.xlsx');
outputDir = fullfile(fileparts(codeDir), 'Analysis_Figures');
addpath(fullfile(codeDir, 'Dependencies', 'ktaub'));

assert(isfile(dataFile), 'Missing input file: %s', dataFile);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

T = readtable(dataFile, 'Sheet', 'Annual area', ...
    'VariableNamingRule', 'preserve');
years = T.Year(:);
assert(isequal(years, (1988:2023)'), ...
    'The input must contain one row per year from 1988 through 2023.');

atlanticArea = T.('Atlantic adjusted (ha)');
pacificArea = T.('Pacific adjusted (ha)');

% Columns are Canada, Atlantic and Pacific. The published national series
% reports adjusted area to the nearest hectare.
area = [round(atlanticArea + pacificArea), atlanticArea, pacificArea];

alpha = 0.05;
windowYears = 10;
lastYears = (min(years) + windowYears - 1:max(years))';
windowMidYears = lastYears - windowYears / 2;

% Calculate the annual-area trend and moving-window acceleration for all
% three geographic series.
for i = 1:size(area, 2)
    result(i).area = area(:, i); %#ok<AGROW>
    [~,~,~,~,~,~,~,result(i).areaSlope,~,~, ...
        result(i).areaLo,result(i).areaHi] = ...
        ktaub([years area(:, i)], alpha, 0);
    intercept = median(area(:, i) - result(i).areaSlope .* years, 'omitnan');
    interceptLo = median(area(:, i) - result(i).areaLo .* years, 'omitnan');
    interceptHi = median(area(:, i) - result(i).areaHi .* years, 'omitnan');
    result(i).areaFit = intercept + result(i).areaSlope .* years;
    result(i).areaFitLo = interceptLo + result(i).areaLo .* years;
    result(i).areaFitHi = interceptHi + result(i).areaHi .* years;
    result(i).areaBandLo = min(result(i).areaFitLo, result(i).areaFitHi);
    result(i).areaBandHi = max(result(i).areaFitLo, result(i).areaFitHi);

    result(i).windowLoss = nan(size(lastYears));
    for j = 1:numel(lastYears)
        use = years >= lastYears(j) - windowYears + 1 & years <= lastYears(j);
        [~,~,~,~,~,~,~,slope] = ktaub([years(use) area(use, i)], alpha, 0);
        result(i).windowLoss(j) = abs(slope);
    end

    [~,~,~,~,~,~,~,result(i).acceleration,~,~, ...
        result(i).accelerationLo,result(i).accelerationHi] = ...
        ktaub([windowMidYears result(i).windowLoss], alpha, 0);
    intercept = median(result(i).windowLoss - ...
        result(i).acceleration .* windowMidYears, 'omitnan');
    interceptLo = median(result(i).windowLoss - ...
        result(i).accelerationLo .* windowMidYears, 'omitnan');
    interceptHi = median(result(i).windowLoss - ...
        result(i).accelerationHi .* windowMidYears, 'omitnan');
    result(i).accelerationFit = intercept + ...
        result(i).acceleration .* windowMidYears;
    result(i).accelerationFitLo = interceptLo + ...
        result(i).accelerationLo .* windowMidYears;
    result(i).accelerationFitHi = interceptHi + ...
        result(i).accelerationHi .* windowMidYears;
    result(i).accelerationBandLo = min( ...
        result(i).accelerationFitLo, result(i).accelerationFitHi);
    result(i).accelerationBandHi = max( ...
        result(i).accelerationFitLo, result(i).accelerationFitHi);
end

%% Create the two figures
fig1 = figure('Color', 'w', 'Position', [60 60 1150 600]);
layout1 = tiledlayout(fig1, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

figED2 = figure('Color', 'w', 'Position', [1 1 1150 1100]);
layoutED2 = tiledlayout(figED2, 2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Each row contains the area and loss-rate axes for Canada, Atlantic and
% Pacific, respectively.
axesPair = {
    nexttile(layout1, 1),   nexttile(layout1, 2);
    nexttile(layoutED2, 1), nexttile(layoutED2, 2);
    nexttile(layoutED2, 3), nexttile(layoutED2, 4)};

areaYLimits = {[85 90], [70.5 73.5], [13.5 16.5]}; % K ha
lossYLimits = {[0 200], [0 120], [30 55]};
panelLetters = {'a', 'b'; 'a', 'b'; 'c', 'd'};
fontSize = 16;
ciColor = [0.90 0.90 0.90];
gridColor = [0.88 0.88 0.88];
gridAlpha = 0.25;

for i = 1:3

    % Annual adjusted area in K ha; trend statistics remain in ha per year.
    ax = axesPair{i, 1};
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    set(ax, 'GridColor', gridColor, 'GridAlpha', gridAlpha);
    ci = patch(ax, [years; flipud(years)], ...
        [result(i).areaBandLo; flipud(result(i).areaBandHi)] / 1000, ciColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.85);
    plot(ax, years, result(i).areaFitLo / 1000, '--', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    plot(ax, years, result(i).areaFitHi / 1000, '--', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    trend = plot(ax, years, result(i).areaFit / 1000, 'k-', 'LineWidth', 1.8);
    points = scatter(ax, years, result(i).area / 1000, 22, 'filled', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    xlabel(ax, 'Year');
    ylabel(ax, 'Tidal marsh area (K ha)');
    xlim(ax, [1988 2023]);
    ylim(ax, areaYLimits{i});
    ax.YRuler.Exponent = 0;
    text(ax, -0.13, 0.99, panelLetters{i, 1}, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', fontSize + 2);
    legend(ax, [points trend ci], ...
        {'Annual area', ...
         sprintf('Trend %.1f ha yr^{-1}', result(i).areaSlope), ...
         sprintf('95%% CI [%.1f, %.1f] ha yr^{-1}', ...
            result(i).areaLo, result(i).areaHi)}, ...
        'Location', 'southwest', 'Box', 'on');

    % Ten-year moving-window loss rate and its Sen trend.
    ax = axesPair{i, 2};
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    set(ax, 'GridColor', gridColor, 'GridAlpha', gridAlpha);
    ci = patch(ax, [windowMidYears; flipud(windowMidYears)], ...
        [result(i).accelerationBandLo; flipud(result(i).accelerationBandHi)], ...
        ciColor, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    plot(ax, windowMidYears, result(i).accelerationFitLo, '--', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    plot(ax, windowMidYears, result(i).accelerationFitHi, '--', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    trend = plot(ax, windowMidYears, result(i).accelerationFit, ...
        'k-', 'LineWidth', 1.8);
    points = scatter(ax, windowMidYears, result(i).windowLoss, 26, 'filled', ...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    xlabel(ax, 'Year');
    ylabel(ax, 'Annual loss rate (ha yr^{-1})');
    xlim(ax, [1988 2023]);
    ylim(ax, lossYLimits{i});
    ax.YRuler.Exponent = 0;
    text(ax, -0.13, 0.99, panelLetters{i, 2}, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', fontSize + 2);
    legend(ax, [points trend ci], ...
        {'Loss rate in 10-yr window', ...
         sprintf('Acceleration %.2f ha yr^{-2}', result(i).acceleration), ...
         sprintf('95%% CI [%.2f, %.2f] ha yr^{-2}', ...
            result(i).accelerationLo, result(i).accelerationHi)}, ...
        'Location', 'southeast', 'Box', 'on');
end

set(findall(fig1, 'Type', 'axes'), 'FontName', 'Arial', ...
    'FontSize', fontSize, 'LineWidth', 1.1, 'Layer', 'top');
set(findall(figED2, 'Type', 'axes'), 'FontName', 'Arial', ...
    'FontSize', fontSize, 'LineWidth', 1.1, 'Layer', 'top');

exportgraphics(fig1, fullfile(outputDir, 'Fig1_TidalMarshTrend.png'), ...
    'Resolution', 600);
exportgraphics(fig1, fullfile(outputDir, 'Fig1_TidalMarshTrend.pdf'), ...
    'ContentType', 'vector');
exportgraphics(figED2, ...
    fullfile(outputDir, 'ExtendedDataFig2_CoastalTrends.png'), ...
    'Resolution', 600);
exportgraphics(figED2, ...
    fullfile(outputDir, 'ExtendedDataFig2_CoastalTrends.pdf'), ...
    'ContentType', 'vector');
