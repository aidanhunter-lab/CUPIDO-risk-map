%% Display location of Antarctic tourist locations on map

%% Preamble
% Include all required (sub)directories within the search path
project = 'CUPIDO-risk-map';
thisFile = which('map_tourism');
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

% Use the IAATO data on Antarctic tourism -- this should be refined later,
% but for now will suffice for a poster plot
filename = fullfile(baseDirectory, 'data','shipping', 'IAATO', ...
    'IAATO-Operator-Site-Visits-by-Activity-2019-2020.csv');
tourism = readtable(filename);

% Exclude any data outside of the map bounding coordinates.
% Map bounding 'vetices' (doesn't really make sense for stereographic projection)
mv = [lon(1), lon(1), lon(2), lon(2), lon(1); ...
    lat(1), lat(2), lat(2), lat(1), lat(1)];
inmap = inpolygon(tourism.Longitude, tourism.Latitude, mv(1,:), mv(2,:)); % find tourist spots within mapped region
tourism = tourism(inmap,:);

% Estimate total visiters at each site -- I think this is more complicated
% than simply adding up the columns, but this will do for a rough approx
% that can be plotted.
fields = fieldnames(tourism);
ind = find(strcmp(fields, 'Aircraft_Flight')):width(tourism);
d = table2array(tourism(:,ind));
visits = sum(d, 2);
tourism.Visits = visits;

% % Identify tourism types
% sampleMethod = unique(plastics.Sampling_Method);
% 
% % Variable to plot
Var = 'Visits';

%% Plot maps

mainTitle = true;
% titleText = {'Tourism locations: 2019-20'};
titleText = {'Tourism locations: 2019-20', 'IAATO site visits'};
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred
titleSize = 13;

% cbarTitleSize = 12; % text sizing for colourbar
% cbarLabelSize = 9;
axisSize = 9;
landColour = .4 .* ones(1,3);
XaxisLocation = 'top';

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

tourism.plotPointSize = nan(height(tourism), 1);
for i = 1:nsizes+1
    if i == 1
        ind = tourism.(Var) <= st(i);
    elseif i == nsizes+1
        ind = st(i-1) < tourism.(Var);
    else
        ind = st(i-1) < tourism.(Var) & tourism.(Var) <= st(i);
    end
    tourism.plotPointSize(ind) = ptSizes(i);
end

% ptSize = 100;

% Choose colours
alpha = 0.65; % transparency
% groups = unique(tourism.(Var));
% ngroups = length(groups);

ptCol = [0.5 0 0.5];

plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
% Create map
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ... 
    'axesLabelSize', axisSize);
hold on

pt = m_scatter(tourism.Longitude, tourism.Latitude, tourism.plotPointSize, ptCol);

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

leg = m_legend_extend([pt5, pt4, pt3, pt2, pt1], '>10000', '>1000', '>100', '>10', '<10');
leg.Title.String = 'visitor numbers';
leg.Title.FontWeight = 'normal';


% manually adjust legend positions

set(pt1, 'Visible', 'off')
set(pt2, 'Visible', 'off')
set(pt3, 'Visible', 'off')
set(pt4, 'Visible', 'off')
set(pt5, 'Visible', 'off')

% save plot
filename = 'tourism_mapPlot_SouthernOcean_IAATO_OperatorSiteVisits2019-20.png';
filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
exportgraphics(plt, filepath)

























