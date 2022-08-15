%% Display Krill-Base data on map

%% Preamble
% Include all required (sub)directories within the search path
project = 'CUPIDO-risk-map';
thisFile = which('map_KrillBase');
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
filename = 'krillbase_data.csv';
krill = readtable(filename);

% Map bounding 'vetices' (doesn't really make sense for stereographic projection)
mv = [lon(1), lon(1), lon(2), lon(2), lon(1); ...
    lat(1), lat(2), lat(2), lat(1), lat(1)];
inmap = inpolygon(krill.LONGITUDE, krill.LATITUDE, mv(1,:), mv(2,:)); % find krill data within mapped region
krill = krill(inmap,:);

% Split date into years, months and days
[yr, mo, day] = datevec(krill.DATE);
krill.year = yr;
krill.month = mo;
krill.day = day;

% Choose a year range -- roughly corresponding to availability of other data
% ymin = 1997;
ymin = min(krill.year);
ymax = max(krill.year);
krill = krill(krill.year >= ymin & krill.year <= ymax,:);

% Use day-time or night-time samples
dayornight = 'both'; % select 'day', 'night', or 'both'
krill = krill(ismember(krill.DAY_NIGHT, {'day', 'night'}),:);
switch dayornight
    case {'day', 'night'}
        krill = krill(strcmp(krill.DAY_NIGHT, dayornight),:);
end

% Variable to plot
% Var = 'NUMBER_OF_KRILL_UNDER_1M2';
Var = 'STANDARDISED_KRILL_UNDER_1M2';

% Omit measuremtns zeros?
omitZeros = true;
switch omitZeros, case true
    krill = krill(krill.(Var) > 0,:);
end


%% Plot maps

% We could plot all data within a year range; all data within separate
% years; or data within months... figure out the best appraoch later.

% For now, let's just plot all data within the selected year range

Months.name = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
Months.nameFull = {'January','February','March','April','May','June','July','August','September','October','November','December'};
Months.index = 1:12;
Months.summerMonths = {'Oct','Nov','Dec','Jan','Feb','Mar'};
Months.winterMonths = {'Apr','May','Jun','Jul','Aug','Sep'};
Months.issummer = ismember(Months.name, Months.summerMonths);
Months.iswinter = ismember(Months.name, Months.winterMonths);

multiPanel = false; % all maps on same figure?
switch multiPanel, case true
    nrows = floor(nmonths ^ 0.5);
    ncols = ceil(nmonths / nrows);
end

mainTitle = true;
titleText = {'Krill trawl surveys'};
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred
titleSize = 13;

% fullTitle = true; % false => abbreviated month
% cbarTitleSize = 12; % text sizing for colourbar
% cbarLabelSize = 9;
axisSize = 9;
landColour = .4 .* ones(1,3);
XaxisLocation = 'top';
% nColourBackTicks = 7; % number of ticks on colourbar
% plotSize = [6 6];
plotSize = [6 7]; % choose same scale as chlorophyll plots so maps are consistent and look good together...
ptScale = 6;
ptSizeMin = 10 .^ (1 / ptScale);
ptSizeMax = 150 .^ (1 / ptScale);

% Choose threshold values to define point sizes
nsizes = 4;
st = 10 .^ (1:nsizes); % threshold values for point sizes
% ptSizes = round(linspace(ptSizeMin, ptSizeMax, nsizes+1));
ptSizes = round(linspace(ptSizeMin, ptSizeMax, nsizes+1) .^ ptScale);


krill.plotPointSize = nan(height(krill), 1);
for i = 1:nsizes+1
    if i == 1
        ind = krill.(Var) <= st(i);
    elseif i == nsizes+1
        ind = st(i-1) < krill.(Var);
    else
        ind = st(i-1) < krill.(Var) & krill.(Var) <= st(i);
    end
    krill.plotPointSize(ind) = ptSizes(i);
end

% Choose colours
alpha = 0.65; % transparency
colDay = [0.8, 0.15, 0];
colNight = [0, 0.45, 0.7];
krill.plotCol = nan(height(krill),3);
ind_day = strcmp(krill.DAY_NIGHT, 'day');
ind_night = strcmp(krill.DAY_NIGHT, 'night');
krill.plotCol(ind_day,:) = repmat(colDay, [sum(ind_day), 1]);
krill.plotCol(ind_night,:) = repmat(colNight, [sum(ind_night), 1]);
ind_day1 = find(ind_day & krill.plotPointSize == ptSizes(3), 1);
ind_night1 = find(ind_night & krill.plotPointSize == ptSizes(3), 1);

plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
% Create map
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ... 
    'axesLabelSize', axisSize);
hold on

pt = m_scatter(krill.LONGITUDE, krill.LATITUDE, krill.plotPointSize, krill.plotCol);

pt_day = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', krill.plotCol(ind_day1,:), ...
    'SizeData', ptSizes(3));
pt_night = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', krill.plotCol(ind_night1,:), ...
    'SizeData', ptSizes(3));

pt1 = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], ...
    'SizeData', ptSizes(1));
pt2 = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], ...
    'SizeData', ptSizes(2));
pt3 = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], ...
    'SizeData', ptSizes(3));
pt4 = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], ...
    'SizeData', ptSizes(4));
pt5 = m_scatter(0, lat(2), 'Marker', 'o', 'MarkerEdgeColor', [0 0 0], ...
    'SizeData', ptSizes(5));

set(pt, {'LineWidth', 'MarkerEdgeAlpha'}, {1, alpha})

switch mainTitle, case true
    m_text(titleLon, titleLat, titleText, ...
        'FontSize', titleSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom')
end

% gc = gca;

leg = m_legend_extend([pt_night, pt_day], 'night', 'day');
leg.Title.String = 'sample time';
leg.Title.FontWeight = 'normal';

leg = m_legend_extend([pt5, pt4, pt3, pt2, pt1], '>10000', '>1000', '>100', '>10', '<10');
leg.Title.String = 'number/m^2';
leg.Title.FontWeight = 'normal';

% manually adjust legend positions

set(pt_day, 'Visible', 'off')
set(pt_night, 'Visible', 'off')
set(pt1, 'Visible', 'off')
set(pt2, 'Visible', 'off')
set(pt3, 'Visible', 'off')
set(pt4, 'Visible', 'off')
set(pt5, 'Visible', 'off')


% save plot
filename = 'KrillBase_mapPlot_SouthernOcean_allSamples.png';
filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
exportgraphics(plt, filepath)






