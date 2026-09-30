%% Figure 5. PCA expansion and tidal-marsh response

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
addpath(fullfile(codeDir, 'Dependencies', 'ktaub'));
rootDir = fileparts(codeDir);
pcaDir = fullfile(rootDir, 'Statistics', 'PCAs-related');
areaDir = fullfile(rootDir, 'Statistics', 'Area-related');
ghmDir = fullfile(rootDir, 'Statistics', 'GHM-related');
outDir = fullfile(rootDir, 'Analysis_Figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end

pcaFile = fullfile(pcaDir, 'Annual_Area_by_PCA_and_Basin_1985_2024.csv');
coastAreaFile = fullfile(areaDir, 'Annual_Area_by_Coast_1988_2023.xlsx');
basinAreaFile = fullfile(ghmDir, 'Annual_Area_by_Basin_1985_2024.csv');

pca = readtable(pcaFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
coastArea = readtable(coastAreaFile, 'VariableNamingRule', 'preserve');
basin = readtable(basinAreaFile, 'VariableNamingRule', 'preserve');

years = (1988:2023)';
nYears = numel(years);
pca.year = double(pca.year);
pca.ValidYear = double(pca.ValidYear);
pca.ID = double(pca.ID);
pca.marsh_area_ha = double(pca.marsh_area_ha);
pca.BIOME = string(pca.BIOME);
pca.Coast = string(pca.Coast);

%% Panel a: annual additions and cumulative coverage
x = [1987; years];
terrestrial = zeros(numel(x), 1);
marine = zeros(numel(x), 1);

terrestrial(1) = sum(pca.marsh_area_ha(pca.year == 1988 & ...
    pca.ValidYear < 1988 & pca.BIOME == "T"), 'omitnan');
marine(1) = sum(pca.marsh_area_ha(pca.year == 1988 & ...
    pca.ValidYear < 1988 & pca.BIOME == "M"), 'omitnan');

for i = 1:nYears
    yr = years(i);
    terrestrial(i + 1) = sum(pca.marsh_area_ha(pca.year == yr & ...
        pca.ValidYear == yr & pca.BIOME == "T"), 'omitnan');
    marine(i + 1) = sum(pca.marsh_area_ha(pca.year == yr & ...
        pca.ValidYear == yr & pca.BIOME == "M"), 'omitnan');
end

coastYears = double(coastArea{:, 1});
totalMarsh = double(coastArea{:, 2}) + double(coastArea{:, 5});
protectedPct = zeros(numel(x), 1);
cumulativeProtected = cumsum(terrestrial + marine);
for i = 1:numel(x)
    yr = max(x(i), 1988);
    denominator = totalMarsh(coastYears == yr);
    protectedPct(i) = 100 * cumulativeProtected(i) / denominator(1);
end

%% Panel b: longest-window basin comparison
basin.ID = double(basin.HYBAS_ID);
basin.year = double(basin.year);
basin.marsh_area_ha = double(basin.marsh_km2) * 100;
basin = basin(ismember(basin.year, years), :);

basinIDs = unique(basin.ID(isfinite(basin.ID)));
[hasBasin, basinRow] = ismember(basin.ID, basinIDs);
[hasYear, basinYear] = ismember(basin.year, years);
basinTotal = accumarray([basinRow(hasBasin & hasYear), basinYear(hasBasin & hasYear)], ...
    basin.marsh_area_ha(hasBasin & hasYear), [numel(basinIDs), nYears], @sum, 0);

pcaTrend = pca(ismember(pca.year, years) & isfinite(pca.ID), :);
[hasPcaBasin, pcaBasinRow] = ismember(pcaTrend.ID, basinIDs);
[hasPcaYear, pcaYearCol] = ismember(pcaTrend.year, years);
pcaByBasin = accumarray([pcaBasinRow(hasPcaBasin & hasPcaYear), ...
    pcaYearCol(hasPcaBasin & hasPcaYear)], ...
    pcaTrend.marsh_area_ha(hasPcaBasin & hasPcaYear), ...
    [numel(basinIDs), nYears], @sum, 0);
outsideByBasin = max(basinTotal - pcaByBasin, 0);

eligible = pcaTrend.ValidYear >= 1995 & pcaTrend.ValidYear <= 2016;
pcaEvents = pcaTrend(eligible, :);
[eventGroup, eventID, eventYear] = findgroups(pcaEvents.ID, pcaEvents.ValidYear);
[~, eventYearCol] = ismember(pcaEvents.year, years);
eventArea = accumarray([eventGroup, eventYearCol], pcaEvents.marsh_area_ha, ...
    [max(eventGroup), nYears], @sum, 0);

preAdvantage = [];
postAdvantage = [];
eventMaxArea = [];
for i = 1:size(eventArea, 1)
    if max(eventArea(i, :)) < 10
        continue
    end

    [hasOutside, outsideRow] = ismember(eventID(i), basinIDs);
    if ~hasOutside || max(outsideByBasin(outsideRow, :)) < 10
        continue
    end

    validYear = eventYear(i);
    window = min(validYear - years(1), years(end) - validYear);
    if window < 7
        continue
    end

    preYears = (validYear - window):(validYear - 1);
    postYears = (validYear + 1):(validYear + window);
    [~, preCols] = ismember(preYears, years);
    [~, postCols] = ismember(postYears, years);
    protectedPre = eventArea(i, preCols);
    protectedPost = eventArea(i, postCols);
    outsidePre = outsideByBasin(outsideRow, preCols);
    outsidePost = outsideByBasin(outsideRow, postCols);

    if max(protectedPre) < 10 || max(protectedPost) < 10 || ...
            max(outsidePre) < 10 || max(outsidePost) < 10
        continue
    end
    if max([protectedPre protectedPost]) == min([protectedPre protectedPost])
        continue
    end

    pcaPreRate = senRate(preYears, protectedPre);
    pcaPostRate = senRate(postYears, protectedPost);
    outsidePreRate = senRate(preYears, outsidePre);
    outsidePostRate = senRate(postYears, outsidePost);

    preAdvantage(end + 1, 1) = pcaPreRate - outsidePreRate; %#ok<SAGROW>
    postAdvantage(end + 1, 1) = pcaPostRate - outsidePostRate; %#ok<SAGROW>
    eventMaxArea(end + 1, 1) = max(eventArea(i, :)); %#ok<SAGROW>
end

%% Plot
fontName = 'Arial';
fontSize = 18;
marineColor = [0 114 178] / 255;
terrColor = [213 94 0] / 255;
pointColor = [0.45 0.45 0.45];
refColor = [0.58 0.35 0.02];

fig = figure('Color', 'w', 'Position', [80 80 1320 570]);

axA = axes(fig, 'Position', [0.072 0.15 0.455 0.73]);
ax = axA; hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
B = bar(ax, x, [terrestrial marine] / 1000, 0.82, 'stacked', 'EdgeColor', 'none');
b1 = B(1); b2 = B(2);
b1.FaceColor = terrColor;
b2.FaceColor = marineColor;
ylabel(ax, 'Annual addition (K ha)');
annualMax = ceil(1.15 * max(terrestrial + marine) / 1000);
ylim(ax, [0 annualMax]);
yticks(ax, 0:1:annualMax);
ax.YColor = [0 0 0];

axRight = axes(fig, 'Position', axA.Position, 'Color', 'none', ...
    'YAxisLocation', 'right', 'XColor', 'none', 'XTick', [], 'YColor', 'k');
hold(axRight, 'on');
plot(axRight, x, protectedPct, '-o', 'Color', [0 0 0], ...
    'LineWidth', 1.7, 'MarkerSize', 4, 'MarkerFaceColor', 'w');
ylabel(axRight, 'Cumulative coverage (%)');
cumulativeMax = ceil(max(protectedPct) / 5) * 5;
ylim(axRight, [0 cumulativeMax]);
yticks(axRight, 0:5:cumulativeMax);

xlim(ax, [1986 2024]);
xticks(ax, [1987 1992 1996 2000 2004 2008 2012 2016 2020]);
xticklabels(ax, {'<1988','1992','1996','2000','2004','2008','2012','2016','2020'});
xtickangle(axA, 0);
xlabel(ax, 'Year');
legendLine = plot(ax, NaN, NaN, '-o', 'Color', 'k', ...
    'LineWidth', 1.7, 'MarkerSize', 4, 'MarkerFaceColor', 'w');
legend(ax, [b1 b2 legendLine], ...
    {'Within terrestrial PCAs', 'Within marine PCAs', 'Cumulative coverage'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', fontSize);

axB = axes(fig, 'Position', [0.663 0.15 0.315 0.73]);
ax = axB; hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
lim = max(1, ceil(max(abs([preAdvantage; postAdvantage]), [], 'omitnan') * 2.2) / 2);
fill(ax, [-lim -lim lim], [-lim lim lim], [0.96 0.92 0.86], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.45);
plot(ax, [-lim lim], [-lim lim], '-', 'Color', refColor, 'LineWidth', 1.3);
xline(ax, 0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
yline(ax, 0, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
sizes = 22 + 190 * sqrt(eventMaxArea ./ max(eventMaxArea, [], 'omitnan'));
scatter(ax, preAdvantage, postAdvantage, sizes, 'filled', ...
    'MarkerFaceColor', pointColor, 'MarkerEdgeColor', 'w', ...
    'MarkerFaceAlpha', 0.62, 'LineWidth', 0.5);
axis(ax, 'square');
xlim(ax, [-lim lim]); ylim(ax, [-lim lim]);
xticks(ax, ceil(-lim):floor(lim)); yticks(ax, ceil(-lim):floor(lim));
xlabel(ax, 'Pre-PCA advantage (% yr^{-1})');
ylabel(ax, 'Post-PCA advantage (% yr^{-1})');

legendAreas = [10 50 500];
h = gobjects(numel(legendAreas), 1);
legendPointColor = 0.62 * pointColor + 0.38;
for i = 1:numel(legendAreas)
    s = 22 + 190 * sqrt(min(legendAreas(i), max(eventMaxArea)) ./ max(eventMaxArea));
    h(i) = plot(ax, NaN, NaN, 'o', 'MarkerSize', sqrt(s), ...
        'MarkerFaceColor', legendPointColor, 'MarkerEdgeColor', 'w');
end
legend(ax, h, compose('%d ha', legendAreas), 'Location', 'northwest', ...
    'Box', 'off', 'FontSize', fontSize);

set(findall(fig, 'Type', 'axes'), 'FontName', fontName, 'FontSize', fontSize, ...
    'LineWidth', 1, 'Layer', 'top', 'XColor', 'k', 'YColor', 'k');
set(findall(fig, 'Type', 'text'), 'FontName', fontName, 'Color', 'k');
set(axRight, 'Position', axA.Position, 'XLim', axA.XLim, ...
    'XColor', 'none', 'XTick', [], 'Box', 'off');
axRight.Toolbar.Visible = 'off';
axesList = [axA axB];
panelLabelOffsets = [0.045 0.058];
for i = 1:2
    ax = axesList(i);
    ax.Toolbar.Visible = 'off';
    annotation(fig, 'textbox', ...
        [ax.Position(1)-panelLabelOffsets(i) ax.Position(2)+ax.Position(4)-0.03 0.03 0.055], ...
        'String', char('a' + i - 1), 'EdgeColor', 'none', 'Margin', 0, ...
        'FontName', fontName, 'FontSize', fontSize + 4, 'FontWeight', 'bold');
end

drawnow;
exportgraphics(fig, fullfile(outDir, 'Fig5_PCAExpansionEffect.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'Fig5_PCAExpansionEffect.pdf'), 'ContentType', 'vector');

function rate = senRate(years, area)
    years = double(years(:));
    logArea = log(max(double(area(:)), 1e-6));
    % Constant series have a zero slope and need no significance test.
    if all(logArea == logArea(1))
        rate = 0;
        return
    end
    [~,~,~,~,~,~,~,slope] = ktaub([years logArea], 0.05, 0);
    rate = 100 * (exp(slope) - 1);
end
