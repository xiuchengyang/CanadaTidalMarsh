%% Extended Data Figure 3. Basin trend distribution
% Plotting section from Analyses/2_Basin/P2_CalculateBasinMarsh.m.

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
figDir = fileparts(codeDir);
srcDir = fullfile(figDir, 'Statistics', 'Area-related');
outDir = fullfile(fileparts(codeDir), 'Analysis_Figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end
addpath(fullfile(codeDir, 'Dependencies', 'hatchfill2_r8'));

T = readtable(fullfile(srcDir, 'Basin_Area_and_Trends_1988_2023.xlsx'), ...
    'VariableNamingRule', 'preserve', 'TextType', 'string');

%% Original 2-ha bins, with a separate no-trend category
trend = T.Trend_ha_per_yr;
binOrder = ["<-10", "-10~-8", "-8~-6", "-6~-4", "-4~-2", "-2~0", ...
    "no trend", "0~2", "2~4", "4~6", "6~8", "8~10", ">10"];
edges = [-inf -10 -8 -6 -4 -2 0 2 4 6 8 10 inf];
numericBins = [1:6 8:13];

binIndex = 7 * ones(height(T), 1);
hasTrend = isfinite(trend) & trend ~= 0;
binIndex(hasTrend) = numericBins(discretize(trend(hasTrend), edges));

category = lower(string(T.category));
stackIndex = 3 * ones(height(T), 1);
stackIndex(startsWith(category, "accelerated")) = 1;
stackIndex(startsWith(category, "decelerated")) = 2;
M = accumarray([binIndex stackIndex], 1, [numel(binOrder) 3]);
assert(sum(M, 'all') == height(T), 'Some basins were not counted.');

%% Colour follows trend magnitude; texture follows acceleration
lossColors = [0.40 0.00 0.00; 0.70 0.00 0.00; 1.00 0.00 0.00; ...
    1.00 0.40 0.40; 1.00 0.60 0.60; 1.00 0.80 0.80];
gainColors = [0.75 0.85 1.00; 0.55 0.70 1.00; 0.25 0.50 1.00; ...
    0.00 0.30 1.00; 0.00 0.20 0.85; 0.00 0.08 0.65];
binColors = [lossColors; 1 1 0.80; gainColors];
fontName = 'Arial';
fontSize = 18;

fig = figure('Color', 'w', 'Position', [80 80 1200 700], 'InvertHardcopy', 'off');
ax = axes(fig, 'Position', [0.09 0.17 0.89 0.80], ...
    'Color', 'w', 'XColor', 'k', 'YColor', 'k');
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 1, ...
    'TickLength', [0 0], 'XMinorTick', 'off', ...
    'GridColor', [0.85 0.85 0.85], 'GridAlpha', 0.15);
ax.Toolbar.Visible = 'off';

bars = bar(ax, M, 'stacked', 'BarWidth', 1);
for k = 1:numel(bars)
    bars(k).FaceColor = 'flat';
    bars(k).CData = binColors;
    bars(k).EdgeColor = [0.6 0.6 0.6];
end
xlim(ax, [0.5 10.5]);
ylim(ax, [0 75]);
xticks(ax, 1:10);
xticklabels(ax, binOrder(1:10));
xtickangle(ax, 10);
yticks(ax, 0:10:70);
xlabel(ax, 'Trend (ha yr^{-1})');
ylabel(ax, 'Number of basins');

%% Cross-hatching for accelerated trends and dots for decelerated trends
hatchLines = hatchfill2(bars(1), 'cross', 'HatchAngle', 45, ...
    'HatchDensity', 28, 'LineWidth', 0.9, 'HatchColor', 'k');

dotDX = 0.12;
dotDY = 0.90;
dotSize = 5;
for i = 1:size(M, 1)
    if M(i, 2) <= 0; continue; end
    xLeft = i - 0.5 + 0.5 * dotDX;
    xRight = i + 0.5 - 0.5 * dotDX;
    yBottom = M(i, 1) + 0.5 * dotDY;
    yTop = M(i, 1) + M(i, 2) - 0.5 * dotDY;
    [xx, yy] = meshgrid(xLeft:dotDX:xRight, yBottom:dotDY:yTop);
    xx(2:2:end, :) = xx(2:2:end, :) + 0.5 * dotDX;
    scatter(ax, xx(:), yy(:), dotSize, 'k', 'filled', ...
        'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
end

%% Pattern legend, matching the manuscript figure
text(ax, 7.95, 72, 'Acceleration', 'FontName', fontName, ...
    'FontSize', fontSize + 1, 'VerticalAlignment', 'middle');
legendPatch = patch(ax, [7.95 8.60 8.60 7.95], [63.5 63.5 68.5 68.5], ...
    'w', 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.9);
hatchfill2(legendPatch, 'cross', 'HatchAngle', 45, ...
    'HatchDensity', 28, 'LineWidth', 0.9, 'HatchColor', 'k');
text(ax, 8.72, 66, 'Accelerated', 'FontName', fontName, ...
    'FontSize', fontSize, 'VerticalAlignment', 'middle');

patch(ax, [7.95 8.60 8.60 7.95], [56.5 56.5 61.5 61.5], ...
    'w', 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.9);
[xx, yy] = meshgrid(8.00:0.065:8.55, 57.0:0.8:61.0);
scatter(ax, xx(:), yy(:), dotSize, 'k', 'filled');
text(ax, 8.72, 59, 'Decelerated', 'FontName', fontName, ...
    'FontSize', fontSize, 'VerticalAlignment', 'middle');

drawnow;
print(fig, fullfile(outDir, 'ExtendedDataFig3_BasinTrendDistribution.png'), '-dpng', '-r300');
exportgraphics(fig, fullfile(outDir, 'ExtendedDataFig3_BasinTrendDistribution.pdf'), ...
    'ContentType', 'vector');
fprintf('Extended Data Fig. 3: %d basins; %d accelerated, %d decelerated, %d without a shift.\n', ...
    sum(M, 'all'), sum(M(:, 1)), sum(M(:, 2)), sum(M(:, 3)));
