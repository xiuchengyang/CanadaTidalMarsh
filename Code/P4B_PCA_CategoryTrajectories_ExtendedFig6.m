%% Extended Data Figure 6. IUCN PCA map and legacy response

clear; clc; close all;

codeDir = fileparts(mfilename('fullpath'));
addpath(fullfile(codeDir, 'Dependencies', 'ktaub'));
rootDir = fileparts(codeDir);
dataDir = fullfile(rootDir, 'Statistics', 'PCAs-related');
outDir = fullfile(rootDir, 'Analysis_Figures');
if ~exist(outDir, 'dir'); mkdir(outDir); end

pcaLongFile = fullfile(dataDir, 'Annual_Area_by_PCA_and_Basin_1985_2024.csv');
basinAreaFile = fullfile(rootDir, 'Statistics', 'GHM-related', ...
    'Annual_Area_by_Basin_1985_2024.csv');
cpcadShp = fullfile(dataDir, 'CPCAD_TidalMarsh_1988_2023', ...
    'CPCAD_TidalMarsh_1988_2023.shp');
basinShp = fullfile(rootDir, 'Statistics', 'Area-related', ...
    'Spatial_boundaries', 'Basin_lev06.shp');

assert(isfile(pcaLongFile), 'Missing input: %s', pcaLongFile);
assert(isfile(basinAreaFile), 'Missing input: %s', basinAreaFile);
assert(isfile(cpcadShp), 'Missing input: %s', cpcadShp);
assert(isfile(basinShp), 'Missing input: %s', basinShp);

groupOrder = ["I_Ia_Ib", "II", "III", "IV"];
groupLabels = ["I (Ia-Ib)", "II", "III", "IV"];
groupShortLabels = ["I", "II", "III", "IV"];
groupColors = [
    0 114 178
    0 158 115
    230 159 0
    204 121 167] ./ 255;
referenceColor = [17 17 17] ./ 255;
focusMaxAge = 25;
fontName = 'Arial';
years = (1988:2023)';

%% Load original overlap data and spatial layers
pcaLong = readtable(pcaLongFile, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
pcaLong.ZONE_ID = regexprep(strtrim(string(pcaLong.ZONE_ID)), '\.0$', '');
pcaLong.Coast = string(pcaLong.Coast);
pcaLong.ID = double(pcaLong.ID);
pcaLong.year = double(pcaLong.year);
pcaLong.ValidYear = double(pcaLong.ValidYear);
pcaLong.marsh_area_ha = double(pcaLong.marsh_area_ha);

fprintf('Reading map shapefiles...\n');
cpcadInfo = shapeinfo(cpcadShp);
cpcad = shaperead(cpcadShp);
basins = shaperead(basinShp);

cpcadZone = regexprep(strtrim(string({cpcad.ZONE_ID}')), '\.0$', '');
catVals = nan(numel(cpcad), 1);
for i = 1:numel(cpcad)
    catVals(i) = str2double(string(cpcad(i).IUCN_CAT));
end

cpcadGroup = strings(numel(cpcad), 1);
cpcadGroup(ismember(round(catVals), [1 2])) = "I_Ia_Ib";
cpcadGroup(round(catVals) == 3) = "II";
cpcadGroup(round(catVals) == 4) = "III";
cpcadGroup(round(catVals) == 5) = "IV";

zoneGroupLookup = unique(table(cpcadZone, cpcadGroup, ...
    'VariableNames', {'ZONE_ID', 'IUCN_Group'}), 'rows');
[hasGroup, groupIdx] = ismember(pcaLong.ZONE_ID, zoneGroupLookup.ZONE_ID);
pcaLong.IUCN_Group = strings(height(pcaLong), 1);
pcaLong.IUCN_Group(hasGroup) = zoneGroupLookup.IUCN_Group(groupIdx(hasGroup));

mapMask = pcaLong.year >= 1988 & pcaLong.year <= 2023 & ...
    pcaLong.marsh_area_ha > 0 & ismember(pcaLong.IUCN_Group, groupOrder);
mapRows = pcaLong(mapMask, :);
marshZones = unique(mapRows.ZONE_ID);

zones = unique(mapRows.ZONE_ID);
zoneCoast = strings(numel(zones), 1);
for i = 1:numel(zones)
    candidates = mapRows.Coast(mapRows.ZONE_ID == zones(i));
    candidates = candidates(strlength(candidates) > 0 & ~ismissing(candidates));
    if ~isempty(candidates)
        u = unique(candidates);
        n = zeros(numel(u), 1);
        for j = 1:numel(u)
            n(j) = sum(candidates == u(j));
        end
        [~, maxIdx] = max(n);
        zoneCoast(i) = u(maxIdx);
    end
end
coastLookup = table(zones, zoneCoast, 'VariableNames', {'ZONE_ID', 'Coast'});

keepCpcad = ismember(cpcadZone, marshZones) & ismember(cpcadGroup, groupOrder);
cpcad = cpcad(keepCpcad);
cpcadZone = cpcadZone(keepCpcad);
cpcadGroup = cpcadGroup(keepCpcad);

[hasCoast, coastIdx] = ismember(cpcadZone, coastLookup.ZONE_ID);
cpcadCoast = strings(numel(cpcad), 1);
cpcadCoast(hasCoast) = coastLookup.Coast(coastIdx(hasCoast));
fprintf('PCA polygons in map: Pacific=%d, Atlantic=%d\n', ...
    sum(cpcadCoast == "Pacific"), sum(cpcadCoast == "Atlantic"));

%% Calculate map counts
countsTbl = table();
for i = 1:numel(groupOrder)
    for coast = ["Pacific", "Atlantic"]
        row = table(coast, groupOrder(i), ...
            sum(cpcadCoast == coast & cpcadGroup == groupOrder(i)), ...
            'VariableNames', {'Coast', 'IUCN_Group', 'N_PCA_Polygons'});
        countsTbl = [countsTbl; row]; %#ok<AGROW>
    end
end

%% Calculate category trajectories
legacyRows = pcaLong(ismember(pcaLong.year, years) & ...
    ismember(pcaLong.IUCN_Group, groupOrder) & isfinite(pcaLong.ValidYear), :);
[legacyGroup, legacyZone, legacyValidYear, legacyIucn] = findgroups( ...
    legacyRows.ZONE_ID, legacyRows.ValidYear, legacyRows.IUCN_Group);
[~, legacyYearCol] = ismember(legacyRows.year, years);
legacyArea = accumarray([legacyGroup, legacyYearCol], legacyRows.marsh_area_ha, ...
    [max(legacyGroup), numel(years)], @sum, 0);

recordGroup = strings(0, 1);
recordAge = zeros(0, 1);
recordEstablishment = zeros(0, 1);
recordCurrent = zeros(0, 1);
recordIndex = zeros(0, 1);
for i = 1:size(legacyArea, 1)
    yr0 = legacyValidYear(i);
    if yr0 < 1988 || yr0 > 2022
        continue
    end
    establishmentArea = legacyArea(i, years == yr0);
    if ~isfinite(establishmentArea) || establishmentArea <= 0
        continue
    end
    maxAge = min(2023 - yr0, focusMaxAge);
    for age = 0:maxAge
        currentArea = legacyArea(i, years == yr0 + age);
        recordGroup(end + 1, 1) = legacyIucn(i); %#ok<SAGROW>
        recordAge(end + 1, 1) = age; %#ok<SAGROW>
        recordEstablishment(end + 1, 1) = establishmentArea; %#ok<SAGROW>
        recordCurrent(end + 1, 1) = currentArea; %#ok<SAGROW>
        recordIndex(end + 1, 1) = currentArea / establishmentArea; %#ok<SAGROW>
    end
end

summaryTbl = table();
for i = 1:numel(groupOrder)
    for age = 0:focusMaxAge
        use = recordGroup == groupOrder(i) & recordAge == age;
        q = weightedQuantiles(recordIndex(use), recordEstablishment(use), ...
            [0.05 0.25 0.75 0.95]);
        areaWeightedIndex = sum(recordCurrent(use)) / sum(recordEstablishment(use));
        row = table(groupOrder(i), age, areaWeightedIndex, ...
            q(1), q(2), q(3), q(4), ...
            'VariableNames', {'IUCN_Group', 'LegacyAgeYears', 'AreaWeightedIndex', ...
            'WeightedP05Index', 'WeightedP25Index', ...
            'WeightedP75Index', 'WeightedP95Index'});
        summaryTbl = [summaryTbl; row]; %#ok<AGROW>
    end
end

area2023 = zeros(numel(groupOrder), 1);
index25 = zeros(numel(groupOrder), 1);
change25 = zeros(numel(groupOrder), 1);
for i = 1:numel(groupOrder)
    activeGroupZones = unique(mapRows.ZONE_ID(mapRows.IUCN_Group == groupOrder(i)));
    use2023 = pcaLong.year == 2023 & pcaLong.IUCN_Group == groupOrder(i) & ...
        ismember(pcaLong.ZONE_ID, activeGroupZones);
    area2023(i) = sum(pcaLong.marsh_area_ha(use2023), 'omitnan');
    row25 = summaryTbl.IUCN_Group == groupOrder(i) & ...
        summaryTbl.LegacyAgeYears == focusMaxAge;
    index25(i) = summaryTbl.AreaWeightedIndex(row25);
    change25(i) = 100 * (index25(i) - 1);
end
ratesTbl = table(groupOrder(:), area2023, index25, change25, ...
    'VariableNames', {'IUCN_Group', 'Area2023_ha', ...
    'IndexAt25Years', 'ChangeAt25Years_pct'});

%% Calculate the marsh-outside-PCA reference
basin = readtable(basinAreaFile, 'VariableNamingRule', 'preserve');
basin.ID = double(basin.HYBAS_ID);
basin.year = double(basin.year);
basin.marsh_area_ha = double(basin.marsh_km2) * 100;
basin = basin(ismember(basin.year, years), :);
basinIDs = unique(basin.ID(isfinite(basin.ID)));
[hasBasin, basinRow] = ismember(basin.ID, basinIDs);
[hasYear, basinYearCol] = ismember(basin.year, years);
basinTotal = accumarray([basinRow(hasBasin & hasYear), basinYearCol(hasBasin & hasYear)], ...
    basin.marsh_area_ha(hasBasin & hasYear), [numel(basinIDs), numel(years)], @sum, 0);

pcaPeriod = pcaLong(ismember(pcaLong.year, years) & isfinite(pcaLong.ID), :);
[hasPcaBasin, pcaBasinRow] = ismember(pcaPeriod.ID, basinIDs);
[hasPcaYear, pcaYearCol] = ismember(pcaPeriod.year, years);
protectedByBasin = accumarray([pcaBasinRow(hasPcaBasin & hasPcaYear), ...
    pcaYearCol(hasPcaBasin & hasPcaYear)], ...
    pcaPeriod.marsh_area_ha(hasPcaBasin & hasPcaYear), ...
    [numel(basinIDs), numel(years)], @sum, 0);
outsideAnnual = sum(max(basinTotal - protectedByBasin, 0), 1)';
validReference = isfinite(years) & isfinite(outsideAnnual) & outsideAnnual > 0;
referenceYears = years(validReference);
referenceLogArea = log(outsideAnnual(validReference));
if numel(referenceYears) < 2
    neverSlope = NaN;
elseif all(referenceLogArea == referenceLogArea(1))
    neverSlope = 0;
else
    [~,~,~,~,~,~,~,neverSlope] = ktaub([referenceYears referenceLogArea], 0.05, 0);
end
neverTbl = table(neverSlope, 'VariableNames', {'sen_log_slope'});

targetCrs = projcrs(3978);
cpcad = reprojectShapes(cpcad, cpcadInfo.CoordinateReferenceSystem, targetCrs, false);
basins = reprojectShapes(basins, [], targetCrs, true);

%% ---------------- Draw combined figure ----------------
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [80 80 1100 1125]);

mapBox = [0.110 0.610 0.780 0.348];
axPacific = axes(fig, 'Position', [0.145 0.616 0.195 0.338]);
drawCoastMapPanel(axPacific, basins, cpcad, cpcadGroup, cpcadCoast, ...
    "Pacific", groupOrder, groupColors);

axAtlantic = axes(fig, 'Position', [0.630 0.616 0.225 0.338]);
drawCoastMapPanel(axAtlantic, basins, cpcad, cpcadGroup, cpcadCoast, ...
    "Atlantic", groupOrder, groupColors);

annotation(fig, 'textbox', [0.078 0.934 0.035 0.030], ...
    'String', 'a', 'EdgeColor', 'none', 'FontName', fontName, ...
    'FontSize', 18.5, 'FontWeight', 'bold', 'Margin', 0);

labelAx = axes(fig, 'Position', mapBox, 'Color', 'none', ...
    'XLim', [0 1], 'YLim', [0 1], 'Visible', 'off');
text(labelAx, 0.025, 0.50, 'Pacific Coast', ...
    'Units', 'normalized', 'FontName', fontName, 'FontSize', 15.0, ...
    'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(labelAx, 0.975, 0.50, 'Atlantic Coast', ...
    'Units', 'normalized', 'FontName', fontName, 'FontSize', 15.0, ...
    'Rotation', 270, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

%% ---------------- Center table in map panel ----------------
tableAx = axes(fig, 'Position', [0.340 0.692 0.285 0.185]);
hold(tableAx, 'on');
axis(tableAx, [0 1 0 1]);
axis(tableAx, 'off');

fontsize = 12.6;
xLine = 0.08;
xCategory = 0.22;
xArea = 0.48;
xPacific = 0.74;
xAtlantic = 0.91;
yHeader = 0.88;
rowYs = [0.68 0.50 0.32 0.14];

text(tableAx, xCategory, yHeader, sprintf('Category\nIUCN'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', fontsize - 0.5);
text(tableAx, xArea, yHeader, sprintf('Area 2023\n(ha)'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', fontsize - 0.8);
text(tableAx, xPacific, yHeader, sprintf('Pacific\nn'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', fontsize - 0.8);
text(tableAx, xAtlantic, yHeader, sprintf('Atlantic\nn'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', fontsize - 0.8);
plot(tableAx, [0.02 0.98], [0.79 0.79], 'Color', [0.60 0.60 0.60], 'LineWidth', 0.6);
for y = [0.59 0.41 0.23]
    plot(tableAx, [0.02 0.98], [y y], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.45);
end

for i = 1:numel(groupOrder)
    group = groupOrder(i);
    rowY = rowYs(i);
    pacificRow = countsTbl.IUCN_Group == group & countsTbl.Coast == "Pacific";
    atlanticRow = countsTbl.IUCN_Group == group & countsTbl.Coast == "Atlantic";
    area2023 = ratesTbl.Area2023_ha(ratesTbl.IUCN_Group == group);

    if any(pacificRow)
        pacificN = countsTbl.N_PCA_Polygons(find(pacificRow, 1, 'first'));
    else
        pacificN = 0;
    end

    if any(atlanticRow)
        atlanticN = countsTbl.N_PCA_Polygons(find(atlanticRow, 1, 'first'));
    else
        atlanticN = 0;
    end

    areaText = char(string(round(area2023(1))));
    areaChunks = {};
    while numel(areaText) > 3
        areaChunks = [{areaText(end-2:end)}, areaChunks]; %#ok<AGROW>
        areaText = areaText(1:end-3);
    end
    areaChunks = [{areaText}, areaChunks];
    areaText = strjoin(areaChunks, ',');

    plot(tableAx, [xLine - 0.05, xLine + 0.05], [rowY rowY], ...
        'Color', groupColors(i, :), 'LineWidth', 3.0, 'Clipping', 'off');
    text(tableAx, xCategory, rowY, groupShortLabels(i), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsize);
    text(tableAx, xArea, rowY, areaText, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsize);
    text(tableAx, xPacific, rowY, sprintf('%d', pacificN), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsize);
    text(tableAx, xAtlantic, rowY, sprintf('%d', atlanticN), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fontsize);
end

annotation(fig, 'rectangle', mapBox, 'Color', [0.47 0.47 0.47], 'LineWidth', 0.85);

%% ---------------- Legacy response panels ----------------
legacyAxes = [
    axes(fig, 'Position', [0.130 0.382 0.350 0.220])
    axes(fig, 'Position', [0.540 0.382 0.350 0.220])
    axes(fig, 'Position', [0.130 0.122 0.350 0.220])
    axes(fig, 'Position', [0.540 0.122 0.350 0.220])
];

neverSlope = neverTbl.sen_log_slope(1);
yValues = [
    summaryTbl.WeightedP05Index
    summaryTbl.WeightedP95Index
    summaryTbl.AreaWeightedIndex
    exp(neverSlope .* (0:focusMaxAge)')];
yValues = yValues(isfinite(yValues));
halfRange = max(abs(yValues - 1.0));
halfRange = max(halfRange * 1.08, 0.05);
yLimits = [1.0 - halfRange, 1.0 + halfRange];

letters = ["b", "c", "d", "e"];
for i = 1:numel(groupOrder)
    ax = legacyAxes(i);
    group = groupOrder(i);
    color = groupColors(i, :);
    data = summaryTbl(summaryTbl.IUCN_Group == group, :);
    data = sortrows(data, 'LegacyAgeYears');
    ratesRow = ratesTbl(ratesTbl.IUCN_Group == group, :);

    x = data.LegacyAgeYears;
    hold(ax, 'on');
    fill(ax, [x; flipud(x)], ...
        [data.WeightedP05Index; flipud(data.WeightedP95Index)], ...
        color, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill(ax, [x; flipud(x)], ...
        [data.WeightedP25Index; flipud(data.WeightedP75Index)], ...
        color, 'FaceAlpha', 0.26, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, x, data.AreaWeightedIndex, 'Color', color, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('IUCN %s', char(groupLabels(i))));

    refX = (0:focusMaxAge)';
    refY = exp(neverSlope .* refX);
    plot(ax, refX, refY, '--', 'Color', referenceColor, 'LineWidth', 1.9, ...
        'DisplayName', 'Marsh outside PCAs');

    xlim(ax, [0 focusMaxAge]);
    ylim(ax, yLimits);
    xticks(ax, 0:5:focusMaxAge);
    grid(ax, 'on');
    box(ax, 'on');
    ax.GridColor = [0.87 0.87 0.87];
    ax.GridAlpha = 0.85;
    ax.LineWidth = 0.8;
    ax.FontSize = 14.5;
    ax.FontName = fontName;
    ax.Layer = 'top';

    if ~ismember(i, [3 4])
        ax.XTickLabel = [];
    else
        xlabel(ax, 'Years after PCA establishment', 'FontSize', 15.6);
    end
    if ismember(i, [1 3])
        ylabel(ax, 'Indexed tidal-marsh area', 'FontSize', 15.6);
    end

    text(ax, -0.075, 1.030, letters(i), 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 18.5, 'FontWeight', 'bold', 'Clipping', 'off');
    text(ax, 0.03, 0.06, ...
        sprintf('25 yr: %.3f (%+.1f%%)', ratesRow.IndexAt25Years(1), ratesRow.ChangeAt25Years_pct(1)), ...
        'Units', 'normalized', 'Color', color, 'FontSize', 14.1, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 14.8);
end

set(findall(fig, '-property', 'FontName'), 'FontName', fontName);
axesList = findall(fig, 'Type', 'axes');
for i = 1:numel(axesList)
    try
        axesList(i).Toolbar.Visible = 'off';
    catch
    end
    try
        disableDefaultInteractivity(axesList(i));
    catch
    end
end

%% ---------------- Export ----------------
pngOut = fullfile(outDir, 'ExtendedDataFig6_PCACategoryTrajectories.png');
pdfOut = fullfile(outDir, 'ExtendedDataFig6_PCACategoryTrajectories.pdf');
drawnow;
figPdf = copyobj(fig, groot);
figPng = copyobj(fig, groot);
set([figPdf, figPng], 'Visible', 'off');

exportgraphics(figPdf, pdfOut, 'ContentType', 'vector');
exportgraphics(figPng, pngOut, 'Resolution', 600);
close([figPdf, figPng]);

fprintf('Wrote:\n  %s\n  %s\n', pngOut, pdfOut);

%% ========================================================================
% Repeated calculations and geometry operations
%% ========================================================================

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

function shapesOut = reprojectShapes(shapesIn, sourceCrs, targetCrs, sourceIsLonLat)
    shapesOut = shapesIn;
    for i = 1:numel(shapesIn)
        valid = isfinite(shapesIn(i).X) & isfinite(shapesIn(i).Y);
        xTarget = nan(size(shapesIn(i).X));
        yTarget = nan(size(shapesIn(i).Y));

        if sourceIsLonLat
            lon = shapesIn(i).X(valid);
            lat = shapesIn(i).Y(valid);
        else
            [lat, lon] = projinv(sourceCrs, shapesIn(i).X(valid), shapesIn(i).Y(valid));
        end

        [xTarget(valid), yTarget(valid)] = projfwd(targetCrs, lat, lon);
        shapesOut(i).X = xTarget;
        shapesOut(i).Y = yTarget;
    end
end

function drawCoastMapPanel(ax, basins, pcas, pcaGroups, pcaCoasts, coast, groupOrder, groupColors)
    hold(ax, 'on');
    for i = 1:numel(basins)
        plot(ax, basins(i).X, basins(i).Y, 'Color', [0.82 0.82 0.82], 'LineWidth', 0.18);
    end

    coastMask = pcaCoasts == coast;
    for groupIdx = 1:numel(groupOrder)
        idx = find(coastMask & pcaGroups == groupOrder(groupIdx));
        for i = idx(:)'
            plot(ax, pcas(i).X, pcas(i).Y, ...
                'Color', groupColors(groupIdx, :), 'LineWidth', 0.62);
        end
    end

    allX = [];
    allY = [];
    coastPcas = pcas(coastMask);
    for i = 1:numel(coastPcas)
        allX = [allX, coastPcas(i).X(isfinite(coastPcas(i).X))]; %#ok<AGROW>
        allY = [allY, coastPcas(i).Y(isfinite(coastPcas(i).Y))]; %#ok<AGROW>
    end

    xmin = min(allX);
    xmax = max(allX);
    ymin = min(allY);
    ymax = max(allY);
    padX = 0.04 * (xmax - xmin);
    padY = 0.04 * (ymax - ymin);

    axis(ax, 'equal');
    xlim(ax, [xmin - padX, xmax + padX]);
    ylim(ax, [ymin - padY, ymax + padY]);
    axis(ax, 'off');
end
