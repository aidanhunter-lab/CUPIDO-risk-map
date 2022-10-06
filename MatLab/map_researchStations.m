%% Display location of research stations on map

%% Preamble
% Include all required (sub)directories within the search path
project = 'CUPIDO-risk-map';
thisFile = which('map_researchStations');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))

%% Map location
coordsTable = readtable('regional bounding coordinates.csv');
disp(coordsTable)

% Choose area from coordsTable -- it's probably useful to vectorise over all areas
area = coordsTable.Location{10};

% Get map coordinates -- these may be specified as input arguments to
% plotBaseMap, and may also be returned as output arguments.
[lon, lat] = plotBaseMap(area, 'createMap', false, 'coordsTable', coordsTable); %, 'coordsTable', coordsTable);


%% Load data

filename = 'COMNAP_Antarctic_Facilities_Master.csv';
filepath = fullfile(baseDirectory, 'data/research stations/COMNAP');
stations = readtable(fullfile(filepath, filename));


%% Create map showing sample locations of various sources

% There are numerous interesting variables in this data table, but those
% most relevant for microplastic release are population size, seasonality
% and station type.
% Use colour to represent population size, and point shape for seasonality.
% Maybe indicate station type with letters.

% Exclude any data outside of the map bounding coordinates.
mv = [lon(1), lon(1), lon(2), lon(2), lon(1); ...
    lat(1), lat(2), lat(2), lat(1), lat(1)];

inmap = inpolygon(stations.Longitude_DD, stations.Latitude_DD, mv(1,:), mv(2,:));
stations = stations(inmap,:);

% Omit stations with zero or NaN population
pop = stations.Peak_Population;
ind = ~isnan(pop) & pop ~= 0;
stations = stations(ind,:);
pop = stations.Peak_Population;

% All countries operating research stations
stations.Operator_primary = strrep(stations.Operator_primary, 'Republic of ', '');
countries = unique(stations.Operator_primary);
ncountries = length(countries);

% Specify a region for all stations
d = stations(cellfun(@(z) isempty(z), stations.Antarctic_Region),:);
% disp(d)
stations.Antarctic_Region(ismember(stations.Record_ID, [62,232]),:) = {'Dronning Maud Land'};
stations.Antarctic_Region(ismember(stations.Record_ID, [13,92]),:) = {'South Orkney Islands'};

% Set point colours according to station population size
lpop = log10(pop);
npopSizes = 5;
sizeClasses = logspace(min(lpop), max(lpop), npopSizes+1);
sizeClasses(2:end-1) = ceil(sizeClasses(2:end-1) / 5) * 5; % round up to nearest 5
% manually adjust to avoid misleading values at the high end of the scale
sizeClasses(3:5) = [30, 100, 200];
groups_size = cell(1, npopSizes);
for i = 1:npopSizes
    groups_size{i} = [num2str(sizeClasses(i)) '-' num2str(sizeClasses(i+1))];
end

fig = figure;
subplot(1,2,1)
histogram(pop)
xlabel('station pop. size')
subplot(1,2,2)
histogram(lpop, log10(sizeClasses))
xlabel('log station pop. size')
close(fig)

% Point colours
cols = cbrewer2('Set1', npopSizes);
stations.plotCol = nan(height(stations), 3);
stations.popSizeClass = cell(height(stations), 1);
for i = 1:height(stations)
    p = stations.Peak_Population(i);
    sc = find(p <= sizeClasses, 1) - 1;
    if p == min(pop), sc = 1; end
    if p == max(pop), sc = npopSizes; end
    stations.plotCol(i,:) = cols(sc,:);
    stations.popSizeClass{i} = groups_size{sc};
end

% Also set point sizes relative to station population size
minSize = 25;
maxSize = 400;
ptSizes = logspace(log10(minSize), log10(maxSize), npopSizes);
% ptSizes = linspace(minSize, maxSize, npopSizes);
stations.plotSize = nan(height(stations), 1);
for i = 1:height(stations)
    p = stations.Peak_Population(i);
    sc = find(p <= sizeClasses, 1) - 1;
    if p == min(pop), sc = 1; end
    if p == max(pop), sc = npopSizes; end
    stations.plotSize(i) = ptSizes(sc);
end

% Set shape according to seasonality
shapes = {'o','s','v','^'};
seasons = unique(stations.Seasonality);
nseasons = length(seasons);
shapes = shapes(1:nseasons);
stations.plotShape = cell(height(stations), 1);
for i = 1:height(stations)
    j = strcmp(stations.Seasonality{i}, seasons);
    stations.plotShape{i} = shapes{j};
end


mainTitle = true;
titleText = {'Research stations'};
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred
titleSize = 13;

% cbarTitleSize = 12; % text sizing for colourbar
% cbarLabelSize = 9;
axisSize = 9;
landColour = .4 .* ones(1,3);
XaxisLocation = 'top';
% nColourBackTicks = 7; % number of ticks on colourbar
% plotSize = [6 6];
plotSize = [6 7]; % use scale that is consistent with other map plots so they look good together on poster
% ptSize = 100;


plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
% Create map
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ... 
    'axesLabelSize', axisSize);
hold on

for i = 1:height(stations)
    m_scatter(stations.Longitude_DD(i), stations.Latitude_DD(i), ...
        'Marker', stations.plotShape(i), 'MarkerEdgeColor', [0,0,0], ...
        'MarkerFaceColor', stations.plotCol(i,:), 'MarkerFaceAlpha', 0.5, ...
        'SizeData', stations.plotSize(i))
end

switch mainTitle, case true
    m_text(titleLon, titleLat, titleText, ...
        'FontSize', titleSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom')
end

% Create legends
% dummy points - colour
pt_col = cell(1, npopSizes);
for i = 1:npopSizes
    d = stations(find(strcmp(stations.popSizeClass, groups_size{i}), 1),:);
    pt_col{npopSizes - i + 1} = m_scatter(0, max(lat), ...
            'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', d.plotCol);
end
groups_col = groups_size;
fgroups_col = flip(groups_col);
% % dummy points - size
% pt_size = cell(1, npopSizes);
% for i = 1:npopSizes
%     d = stations(find(stations.plotSize == ptSizes(i), 1),:);
%     pt_size{npopSizes - i + 1} = m_scatter(0, max(lat), ...
%             'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', [1,1,1], ...
%             'SizeData', d.plotSize);
% end
% groups_size = cell(1, npopSizes);
% for i = 1:npopSizes
%     groups_size{i} = [num2str(sizeClasses(i)) '-' num2str(sizeClasses(i+1))];
% end
% fgroups_size = flip(groups_size);
% dummy points - size
pt_sh = cell(1, nseasons);
for i = 1:nseasons
    d = stations(find(strcmp(stations.Seasonality, seasons{i}), 1),:);
    pt_sh{nseasons - i + 1} = m_scatter(0, max(lat), ...
            'Marker', d.plotShape, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', [1,1,1]);
end
groups_sh = seasons;
fgroups_sh = flip(groups_sh);
% legends
leg_col = m_legend_extend([pt_col{:}], fgroups_col{:});
legPos_col = get(leg_col, 'Position');
set(leg_col, 'Position', [0.1, 0.1, legPos_col(3), legPos_col(4)]);
legPos_col = get(leg_col, 'Position');
cellfun(@(z) set(z, 'Visible', 'off'), pt_col) % remove dummy points from plot
set(leg_col.Title, 'String', 'population size')
leg_col.Position(3) = 1.1 * leg_col.Position(3);

leg_sh = m_legend_extend([pt_sh{:}], fgroups_sh{:});
legPos_sh = get(leg_sh, 'Position');
set(leg_sh, 'Position', [legPos_col(1) + legPos_col(3) + 0.2 * legPos_col(3), 0.1, legPos_sh(3), legPos_sh(4)]);
cellfun(@(z) set(z, 'Visible', 'off'), pt_sh) % remove dummy points from plot
leg_sh.Position(3) = 1.13 * leg_sh.Position(3);

% save plot
filename = 'researchStation_mapPlot.png';
filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
exportgraphics(plt, filepath)

