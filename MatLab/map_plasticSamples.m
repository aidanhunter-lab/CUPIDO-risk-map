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

plastics = load_plastic_conc_data;

%% Filter data

[abundance, sources, nsources] = filter_plastic_conc_data(plastics, 'saveData', true);

% For each data source, the abundance field should contain the information
% that we want to plot here -- longitude-latitude, abundance, and data 
% source. A nice plot will show abundance as point size and source as point
% colour or shape. 
% The initial basic plot will simply display all samples and a measure of
% total microplastic abundance, without disciminating between plastic
% categories, polymers, colours or any other plastic property. However, the
% plot will exclude macroplastic samples as these are not bioavailble to
% zooplankton.




%% Summary boxplots
cols = cbrewer2('Set2', nsources); % Choose colour for each data source
Cols = table(sources, cols);
Cols.Properties.VariableNames = {'Source', 'Col'};

alpha = 0.5; % transparency
logScale = true;
Depths = {'all', 'surface','subsurface'};
Depth = Depths{2};
Types = {'total', 'fragment', 'fibre'};
Type = Types{3};

for p = 1:length(Types) % separate plot for each plastic type
    Type = Types{p};

    d = abundance;
    switch logScale, case true, d.Value(d.Value == 0) = nan; end
    j = strcmp(d.Type, Type);
    d = d(j,:);
    switch Depth
        case 'all', j = true(height(d),1);
        case 'surface', j = ismember(d.Depth, {'surface', '5m', '<1m'});
        case 'subsurface', j = ismember(d.Depth, {'subsurface'});
    end
    d = d(j,:);
    Variables = unique(d.Variable);
    nVariables = length(Variables);
    useUnit = unique(d(:,{'Variable', 'Unit'}));

    pltName = ['plt_' Type];
    assignin('base', pltName, figure)
    set(evalin('base', pltName), {'Units', 'Position'}, {'inches', [0 0 6 * nVariables, 6]})

    for i = 1:nVariables
        subplot(1, nVariables, i)
        d_ = d(strcmp(d.Variable, Variables{i}),:);
        x = d_.Value;
        g = categorical(d_.Source);
        b = boxchart(x, 'GroupByColor', g);
        bs = arrayfun(@(z) z.DisplayName, b, 'UniformOutput', false);
        for j = 1:length(bs)
            k = Cols.Col(strcmp(Cols.Source, bs{j}),:);
            set(b(j), {'BoxFaceColor', 'BoxLineColor', 'MarkerColor', 'BoxFaceAlpha'}, ...
                {k, k, k, alpha})
        end
        set(gca, {'XTick', 'YScale'}, {[], 'log'})
        pu = useUnit.Unit{strcmp(useUnit.Variable, Variables{i})};

        if ~isempty(str2double(pu(end))), pu = [pu(1:end-1) '^' pu(end)]; end

        switch Variables{i}
            case 'concentration', ylab = ['concentration (', pu, ')'];
            case 'density', ylab = ['density (', pu, ')'];
            case 'massDensity', ylab = ['mass density (', pu, ')'];
        end
        ylabel(ylab)
        legend('Location', 'northwest');
    end
    switch Type
        case 'total', sgtitle('Total microplastic: in situ samples')
        otherwise, sgtitle(['Microplastic ' Type 's: in situ samples'])
    end

end





%% Create map showing sample locations of various sources

% coordsTable = readtable('regional bounding coordinates.csv');
% disp(coordsTable)
% 
% % Choose area from coordsTable -- it's probably useful to vectorise over all areas
% area = coordsTable.Location{10};
% 
% % Get map coordinates -- these may be specified as input arguments to
% % plotBaseMap, and may also be returned as output arguments.
% [lon, lat] = plotBaseMap(area, 'createMap', false, 'coordsTable', coordsTable); %, 'coordsTable', coordsTable);

% Exclude any data outside of the map bounding coordinates.
mv = [lon(1), lon(1), lon(2), lon(2), lon(1); ...
    lat(1), lat(2), lat(2), lat(1), lat(1)];

inmap = inpolygon(abundance.Longitude, abundance.Latitude, mv(1,:), mv(2,:));
abundance = abundance(inmap,:);


mainTitle = true;
titleText = {'Microplastic samples'};
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

% % Choose colours
groups = unique(abundance.Source, 'stable');
ngroups = length(groups);
% cols = cbrewer2('Set1', ngroups);
% %     'Accent',   'qual'; ...
% %     'Dark2',    'qual'; ...
% %     'Paired',   'qual'; ...
% %     'Pastel1',  'qual'; ...
% %     'Pastel2',  'qual'; ...
% %     'Set1',     'qual'; ...
% %     'Set2',     'qual'; ...
% %     'Set3',     'qual'; ...

% alpha = 0.5; % transparency

abundance.plotCol = nan(height(abundance), 3);
for i = 1:ngroups
    j = strcmp(abundance.Source, groups{i});
    k = Cols.Col(strcmp(Cols.Source, groups{i}),:);
    abundance.plotCol(j,:) = repmat(k, sum(j), 1);
end

plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
% Create map
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ... 
    'axesLabelSize', axisSize);
hold on

for i = 1:height(abundance)
    m_scatter(abundance.Longitude(i), abundance.Latitude(i), ...
        'MarkerEdgeColor', abundance.plotCol(i,:), ...
        'MarkerFaceColor', abundance.plotCol(i,:), 'MarkerFaceAlpha', alpha)
end

switch mainTitle, case true
    m_text(titleLon, titleLat, titleText, ...
        'FontSize', titleSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom')
end

% Create legend
pt = cell(1, ngroups); % dummy points for each source
for i = 1:ngroups
    pt{ngroups - i + 1} = m_scatter(0, lat(2), 'MarkerEdgeColor', cols(i,:), ... 
        'MarkerFaceColor', cols(i,:));
end
fgroups = flip(groups);
leg = m_legend_extend([pt{:}], fgroups{:});
legPos = get(leg, 'Position');
set(leg, 'Position', [0.1, 0.1, legPos(3), legPos(4)]);
cellfun(@(z) set(z, 'Visible', 'off'), pt) % remove dummy points from plot

% save plot
filename = 'plastic_mapPlot_SouthernOcean_colouredBySource.png';
filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
exportgraphics(plt, filepath)

