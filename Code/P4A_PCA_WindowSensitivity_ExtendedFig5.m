%% Extended Data Figure 5. PCA window-length sensitivity

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
addpath(fullfile(codeDir, 'Dependencies', 'ktaub'));
rootDir = fileparts(codeDir);
dataDir = fullfile(rootDir, 'Statistics', 'PCAs-related');
outDir = fullfile(rootDir, 'Analysis_Figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end

pcaFile = fullfile(dataDir, 'Annual_Area_by_PCA_and_Basin_1985_2024.csv');
pca = readtable(pcaFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');

years = (1988:2023)';
windowLengths = (4:15)';
pca.ZONE_ID = regexprep(strtrim(string(pca.ZONE_ID)), '\.0$', '');
pca.ValidYear = double(pca.ValidYear);
pca.year = double(pca.year);
pca.marsh_area_ha = double(pca.marsh_area_ha);
pca = pca(ismember(pca.year, years) & isfinite(pca.ValidYear), :);

[zoneGroup, zoneID, validYear] = findgroups(pca.ZONE_ID, pca.ValidYear);
[~, yearCol] = ismember(pca.year, years);
zoneArea = accumarray([zoneGroup, yearCol], pca.marsh_area_ha, ...
    [max(zoneGroup), numel(years)], @sum, 0);

rates = [];
weights = [];
for i = 1:size(zoneArea, 1)
    yr0 = validYear(i);
    if yr0 < 1988 || yr0 + max(windowLengths) > 2023
        continue
    end

    establishmentArea = zoneArea(i, years == yr0);
    if ~isfinite(establishmentArea) || establishmentArea < 10
        continue
    end

    currentRates = nan(1, numel(windowLengths));
    for j = 1:numel(windowLengths)
        postYears = (yr0 + 1):(yr0 + windowLengths(j));
        postArea = zoneArea(i, ismember(years, postYears));
        logArea = log(max(double(postArea(:)), 1e-6));
        if all(logArea == logArea(1))
            slope = 0;
        else
            [~,~,~,~,~,~,~,slope] = ktaub([postYears(:) logArea], 0.05, 0);
        end
        currentRates(j) = 100 * (exp(slope) - 1);
    end

    if all(isfinite(currentRates))
        rates(end + 1, :) = currentRates; %#ok<SAGROW>
        weights(end + 1, 1) = establishmentArea; %#ok<SAGROW>
    end
end

meanRate = nan(numel(windowLengths), 1);
p10 = meanRate; p25 = meanRate; p75 = meanRate; p90 = meanRate;
for j = 1:numel(windowLengths)
    meanRate(j) = sum(rates(:, j) .* weights) / sum(weights);
    q = weightedQuantiles(rates(:, j), weights, [0.10 0.25 0.75 0.90]);
    p10(j) = q(1); p25(j) = q(2); p75(j) = q(3); p90(j) = q(4);
end

fontName = 'Arial';
fontSize = 13;
bandColor = [110 181 168] / 255;

fig = figure('Color', 'w', 'Position', [80 80 610 460]);
ax = axes(fig); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
fill(ax, [windowLengths; flipud(windowLengths)], [p10; flipud(p90)], bandColor, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.16);
fill(ax, [windowLengths; flipud(windowLengths)], [p25; flipud(p75)], bandColor, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.34);
plot(ax, windowLengths, meanRate, '-', 'Color', [0 0 0], 'LineWidth', 2.1);
yline(ax, 0, '--', 'Color', [0.43 0.43 0.43], 'LineWidth', 1);

xlim(ax, [min(windowLengths) max(windowLengths)]);
ylim(ax, [-0.5 0.5]);
xticks(ax, windowLengths);
xlabel(ax, 'Post-PCA establishment window (yr)');
ylabel(ax, 'Tidal-marsh change rate (% yr^{-1})');
legend(ax, {'Weighted 10-90%', 'Weighted 25-75%', 'Weighted mean'}, ...
    'Location', 'northeast', 'Box', 'off');

set(ax, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', 1, ...
    'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
set(findall(fig, 'Type', 'text'), 'FontName', fontName, 'Color', 'k');

exportgraphics(fig, fullfile(outDir, 'ExtendedDataFig5_PCAWindowSensitivity.png'), ...
    'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'ExtendedDataFig5_PCAWindowSensitivity.pdf'), ...
    'ContentType', 'vector');

function q = weightedQuantiles(values, weights, probabilities)
    valid = isfinite(values) & isfinite(weights) & weights > 0;
    values = values(valid);
    weights = weights(valid);
    [values, order] = sort(values);
    weights = weights(order);
    cumulative = cumsum(weights);
    targets = probabilities .* cumulative(end);
    targets = min(max(targets, cumulative(1)), cumulative(end));
    q = interp1(cumulative, values, targets, 'linear');
end
