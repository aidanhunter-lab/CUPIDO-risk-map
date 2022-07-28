%% Display location of plastic samples on map

%% Preamble
% Include all required (sub)directories within the search path
project = 'CUPIDO-risk-map';
thisFile = which('map_plasticSamples');
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

% Waller (2017) providews a list of plastic samples with lat-lon coords.
% Just use this for now -- I can include more later...
filename = fullfile(baseDirectory, 'data','plastic_quantity', 'Waller_2017', 'tableS1.csv');
plastics = readtable(filename);

% Exclude any data outside of the map bounding coordinates.
% Map bounding 'vetices' (doesn't really make sense for stereographic projection)
mv = [lon(1), lon(1), lon(2), lon(2), lon(1); ...
    lat(1), lat(2), lat(2), lat(1), lat(1)];
inmap = inpolygon(plastics.Long, plastics.Lat, mv(1,:), mv(2,:)); % find plastic data within mapped region
plastics = plastics(inmap,:);

% Identify plastic types -- macro/micro
plasticTypes = unique(plastics.Type);

% Identify sampling method
sampleMethod = unique(plastics.Sampling_Method);

% Variable to plot
Var = 'Type';

%% Plot maps

mainTitle = true;
titleText = {'Plastic samples'};
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred
titleSize = 13;

% cbarTitleSize = 12; % text sizing for colourbar
% cbarLabelSize = 9;
axisSize = 9;
landColour = .4 .* ones(1,3);
XaxisLocation = 'top';
% nColourBackTicks = 7; % number of ticks on colourbar
plotSize = [6 7]; % include a little extra height if using multi-line titles
% ptSizeMin = 10;
% ptSizeMax = 200;
ptSize = 100;

% Choose colours
alpha = 0.5; % transparency
groups = unique(plastics.(Var));
ngroups = length(groups);

colMacro = [0.4 0.8 0];
colMicro = [0.8 0 0.4];

plastics.plotCol = nan(height(plastics),3);
ind_macro = strcmp(plastics.(Var), 'Macro');
ind_micro = strcmp(plastics.(Var), 'Micro');
plastics.plotCol(ind_macro,:) = repmat(colMacro, [sum(ind_macro), 1]);
plastics.plotCol(ind_micro,:) = repmat(colMicro, [sum(ind_micro), 1]);
ind_macro1 = find(ind_macro, 1);
ind_micro1 = find(ind_micro, 1);

plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
% Create map
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ... 
    'axesLabelSize', axisSize);
hold on

for i = height(plastics):-1:1
    m_scatter(plastics.Long(i), plastics.Lat(i), 'SizeData', ptSize, ...
        'MarkerEdgeColor', plastics.plotCol(i,:), 'MarkerFaceColor', plastics.plotCol(i,:), ...
        'MarkerFaceAlpha', alpha);
end

pt_macro = m_scatter(0, lat(2), 'SizeData', ptSize, ...
        'MarkerEdgeColor', plastics.plotCol(ind_macro1,:), 'MarkerFaceColor', plastics.plotCol(ind_macro1,:), ...
        'MarkerFaceAlpha', alpha);
pt_micro = m_scatter(0, lat(2), 'SizeData', ptSize, ...
        'MarkerEdgeColor', plastics.plotCol(ind_micro1,:), 'MarkerFaceColor', plastics.plotCol(ind_micro1,:), ...
        'MarkerFaceAlpha', alpha);

switch mainTitle, case true
    m_text(titleLon, titleLat, titleText, ...
        'FontSize', titleSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom')
end

% gc = gca;

% The m_legend.m function from package m_map has been modified -- see the utility directory
leg = m_legend_extend([pt_macro, pt_micro], 'macro', 'micro');
leg.Title.String = 'plastic size';
leg.Title.FontWeight = 'normal';

% manually adjust legend positions

set(pt_macro, 'Visible', 'off')
set(pt_micro, 'Visible', 'off')


% save plot
filename = 'plastic_mapPlot_SouthernOcean_Waller2017.png';
filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
exportgraphics(plt, filepath)

























