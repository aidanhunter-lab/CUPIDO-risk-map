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
sources = fieldnames(plastics);
nsources = length(sources);

% For each data source, the abundance field should contain the information
% that we want to plot here -- longitude-latitude, abundance, and data 
% source. A nice plot will show abundance as point size and source as point
% colour or shape. 
% The initial basic plot will simply display all samples and a measure of
% total microplastic abundance, without disciminating between plastic
% categories, polymers, colours or any other plastic property. However, the
% plot will exclude macroplastic samples as these are not bioavailble to
% zooplankton.


% What measurement units were used by each source?
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    units.(source) = unique(d.Unit, 'stable')';
end
disp(units)
% Volumetric concentration 'pieces/m^3' is most common unit; some data use
% area density 'pieces/km^2'
varOpts = {'concentration', 'density', 'massDensity'}; %, 'all'}; % filtering options
Units.concentration = {'pieces/m3', 'items/m3', 'particles/m3', 'number/m3'};
Units.density = {'pieces/km2', 'items/km2', 'number/km2'};
Units.massDensity = {'g/km2'};
% Units.all = struct2table(units);
% Units.all = unique(Units.all{1,:});
% Units.all = Units.all(~ismember(Units.all, ...
%     {char(), 'g', 'm2' , 'm3', 'number', 'pieces'}));
Units_ = struct2cell(Units);
useUnit = cellfun(@(z) z(1), Units_);

% Include a Variable column in each data set to indicate measurement type,
% omitting rows with units not listed in Units
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    for j = height(d):-1:1
        k = cellfun(@(z) ismember(d.Unit{j}, z), Units_);
        if any(k)
            d.Variable(j) = varOpts(k);
            d.Unit(j) = useUnit(k);
        else
            d.Variable(j) = {''};
        end
    end
    d = d(~strcmp(d.Variable, ''),:);
    plastics.(source).abundance = d;
end



% Regularise the data and combine into a single table
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    fields = d.Properties.VariableNames;
    % Regularise the Station variable
    j = ismember(fields, 'station');
    if any(j), fields{j} = 'Station'; end
    j = ismember('Station', fields);
    if ~j
        k = contains(fields, 'Sample');
        if any(k)
            if sum(k) == 1
                fields{k} = 'Station';
            else
                warning(['Multiple columns contain keyword "Sample" in ' source ' data. Rename columns before running loop.'])
            end
        end
    end
    d.Properties.VariableNames = fields;
    j = ismember('Station', fields);
    if ~j, error(['No variable has been (re)named as "Station" in ' source ' data!']); end
    % Code unique stations as numeric integer lists: 1,2,3,...
    x = d(:,'Station');
    Station = unique(x.Station, 'stable');
    nx = length(Station);
    StationNew = (1:nx)';
    ux = table(Station, StationNew);
    x = join(x, ux);
    d.Station = x.StationNew;
    d = movevars(d, 'Station', 'Before', fields{1});
    fields = d.Properties.VariableNames;
    % Some data sources specify a station Location
    j = ismember(fields, {'Location', 'location', 'Region', 'region'});
    if ~any(j)
        d.Location = cell(height(d), 1);
    else
        d.Properties.VariableNames{j} = 'Location';
    end
    d = movevars(d, 'Location', 'After', 'Station');
    fields = d.Properties.VariableNames;
    % Date
    j = ismember(fields, {'Date', 'date'});
    if ~any(j)
        d.Date = NaT(height(d), 1);
%         d.Date = cell(height(d), 1);
    else
        d.Properties.VariableNames{j} = 'Date';
    end
    d = movevars(d, 'Date', 'After', 'Location');
    fields = d.Properties.VariableNames;
    % Regularise Lon_lat coordinates
    % First Longitude
    j = ismember('Longitude', fields);
    if ~j
        k = ismember(fields, {'longitude', 'Lon', 'lon'});
        if any(k)
            fields(k) = 'Longitude';
        else
            error(['There is no "Longitude" variable in ' source ' data! Check variable names.'])
        end
    end
    d.Properties.VariableNames = fields;
    if iscell(d.Longitude), d.Longitude = cell2mat(d.Longitude); end
    d = movevars(d, 'Longitude', 'After', 'Date');
    fields = d.Properties.VariableNames;
    % Now Latitude    
    j = ismember('Latitude', fields);
    if ~j
        k = ismember(fields, {'latitude', 'Lat', 'lat'});
        if any(k)
            fields(k) = 'Latitude';
        else
            error(['There is no "Latitude" variable in ' source ' data! Check variable names.'])
        end
    end
    d.Properties.VariableNames = fields;
    if iscell(d.Latitude), d.Latitude = cell2mat(d.Latitude); end
    d = movevars(d, 'Latitude', 'After', 'Longitude');
    fields = d.Properties.VariableNames;
    % Regularise Depth
    j = ismember(fields, {'Depth', 'depth', 'Dep', 'dep'});
    if ~any(j)
        d.Depth = cell(height(d), 1);
    else
        d.Properties.VariableNames{j} = 'Depth';
    end
    d = movevars(d, 'Depth', 'After', 'Latitude');
%     fields = d.Properties.VariableNames;
    % Regularise Variable -- this is already selected above
%     Variable = Variables{cellfun(@(z) ismember(d.Unit(1), z), useUnits_)};
%     j = ismember(fields, {'Variable', 'variable'});
%     if ~any(j)
%         d.Variable = repmat({Variable}, height(d), 1);
%     else
%         d.Properties.VariableNames{j} = 'Variable';
%         d.Variable = repmat({Variable}, height(d), 1);
%     end
    d = movevars(d, 'Variable', 'After', 'Depth');
    % Unit was filtered above
    d = movevars(d, 'Unit', 'After', 'Variable');
    fields = d.Properties.VariableNames;
    % Measurement Statistic
    j = ismember(fields, {'Stat', 'stat', 'Statistic', 'statistic', 'Measure', 'measure'});
    if ~any(j)
        d.Stat = repmat({'mean'}, height(d), 1);
    else
        d.Properties.VariableNames{j} = 'Stat';
    end
    x = d.Stat;
    j = ismember(x, {'Min', 'min', 'Minimum', 'minimum'});
    x(j) = {'min'};
    j = ismember(x, {'Max', 'max', 'Maximum', 'maximum'});
    x(j) = {'max'};
    j = ismember(x, {'Mean', 'mean', 'Average', 'average'});
    x(j) = {'mean'};
    j = ismember(x, {'SD', 'sd', 'StDev', 'stdev'});
    x(j) = {'sd'};
    d.Stat = x;
    d = movevars(d, 'Stat', 'After', 'Unit');
    fields = d.Properties.VariableNames;
    % Value
    j = ismember(fields, {'Value', 'value'});
    if ~any(j)
        warning(['There may not be a "Value" column in ' source ' data! Check the variable names.'])
        d.Variable = cell(height(d), 1);
    else
        d.Properties.VariableNames{j} = 'Value';
    end
    d = movevars(d, 'Value', 'After', 'Stat');
    plastics.(source).abundance = d;
end

% Only use microplastic measures
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    fields = d.Properties.VariableNames;
    j = ismember(fields, {'Size','size'});
    if ~any(j), continue; end
    j = ismember(d{:,j}, {'Micro', 'micro'});
    plastics.(source).abundance = d(j,:);
end

% Include a "Source" column in each abundance data set
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    source_ = source;
    j = strfind(source, '_');
    if ~isempty(j)
        source_(j) = ' ';
        source_ = [source_(1:end-4) '(' source_(end-3:end) ')'];
    end
    d.Source = repmat({source_}, height(d), 1);
    d = movevars(d, 'Source', 'Before', 'Station');
    plastics.(source).abundance = d;
end

% Some data sources do not report the measurement means, but only the min 
% and max measurements at each station. This is usually because the data
% were reported as within some range [min, max]. Let's estimate mean values
% for those data where only min/max is reported simply by taking midpoints.

% useStat = 'max';
useStat = 'mean';
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    stats = unique(d.Stat, 'stable');
    useStatExists = ismember(useStat, stats);
    if useStatExists
        d = d(strcmp(d.Stat, useStat),:);
        plastics.(source).abundance = d;
        continue
    else
        switch useStat
            case 'mean'
                j = ismember({'min','max'}, stats);
                if all(j)
                    stations = unique(d.Station);
                    nstations = length(stations);
                    inc = height(d) / nstations;
                    d_ = d(inc:inc:height(d),:);
                    d_.Stat = repmat({useStat}, height(d_), 1);
                    for k = 1:length(stations)
                        x = d.Value(d.Station == k & ismember(d.Stat, {'min','max'}));
                        if x(1) > 0
                            x = exp(0.5 * sum(log(x)));
                        elseif x(1) == 0
                            x = 0.5 * x(2);
                        end
                        d_.Value(d_.Station ==k) = x;
                    end
                    d = d_;
                else
                    warning(['Data from ' source ' has been omitted because useStat is unmatched and the Stat column does not contain rows labelled "min" and "max".'])
                end
            case 'max'
                j = strcmp(d.Stat, 'mean');
                if ~any(j)
                    warning(['Data from ' source ' has been omitted because useStat is unmatched and the Stat column does not contain rows labelled "mean".'])
                end
                d = d(j,:);
        end
        plastics.(source).abundance = d;
    end
end % The abundance data should now contain a single measurement per station

% Include info on plastic type (particle, fibre, etc).
for i = 1:nsources
    % Combine info on plastic category into the abundance data
    source = sources{i};
    abundance = plastics.(source).abundance;
    categories = plastics.(source).categories;
    fields = categories.Properties.VariableNames;
    Stats = {'Measure', 'Stat'}; % We're only interested in mean values for category
    Stat = ismember(Stats, fields);
    if any(Stat)
        y = categories.(Stats{Stat});
        j = ismember(y, {'mean', 'Mean', 'average', 'Average'});
        categories = categories(j,:);
    end
    Types = unique(categories.Type, 'stable');
    nTypes = length(Types);
    abundance.Type = repmat({'total'}, height(abundance), 1);
    abundance = movevars(abundance, 'Type', 'Before', 'Variable');
    Depths = unique(abundance.Depth, 'stable');
    nDepths = length(Depths);
    accountForDepth = ismember('Depth', categories.Properties.VariableNames);
    if nDepths == 1, accountForDepth = false; end
    if accountForDepth
        cDepths = unique(categories.Depth, 'stable');
        if length(cDepths) ~= nDepths || ~all(strcmp(cDepths, Depths))
            warning("There's a mismatch in the depths stored in plastic abundance and plastic category data.")
        end
    end
    clearvars x
    x.total = abundance;
    for j = 1:nTypes
        Type = Types{j};
        x.(Type) = abundance;
        x.(Type).Type = repmat({Type}, height(abundance), 1);
        if accountForDepth
            for k = 1:nDepths
                Depth = Depths{k};
                y = categories(strcmp(categories.Type, Type) & ...
                    strcmp(categories.Depth, Depth),:);
                ind = strcmp(x.(Type).Depth, Depth);
                x.(Type).Value(ind) = y.Value .* x.(Type).Value(ind);
            end
        else
            x_ = x.(Type)(:,{'Station', 'Unit', 'Value'});
            y = categories(strcmp(categories.Type, Type),:);
            if ~ismember('Station', fieldnames(y))
                x.(Type).Value = y.Value .* x.(Type).Value;
            else
                y_ = y(:,{'Station', 'Value'});
                x_ = x.(Type)(:,{'Station', 'Unit', 'Value'});
                x_ = innerjoin(x_, y_, 'Keys', 'Station');
                x.(Type).Value = x_.Value_x_ .* x_.Value_y_;
            end
%             y = categories(strcmp(categories.Type, Type),:);
%             x.(Type).Value = y.Value .* x.(Type).Value;
        end
        if ismember(y.Unit, {'percent','Percent'})
            x.(Type).Value = 0.01 .* x.(Type).Value;
        end
    end
    abundance = struct2cell(x);
    plastics.(source).abundance = vertcat(abundance{:});
end

% Regularise the naming convention for plastic types
for i = 1:nsources
    source = sources{i};
    abundance = plastics.(source).abundance;
    Type = abundance.Type;
    Fibres = {'Line','line','Fibre','fibre','microfibre',...
        'Microfibre','microFibre','MicroFibre'};
    isFibre = ismember(Type, Fibres);
    Type(isFibre) = repmat({'fibre'}, sum(isFibre), 1);
    Fragments = {'fragment','Fragment','flake','Flake','granule',...
        'Granule','sphere','Sphere','microplastic','Microplastic',...
        'microPlastic','MicroPlastic', 'particle', 'Particle'};
    isFragment = ismember(Type, Fragments);
    Type(isFragment) = repmat({'fragment'}, sum(isFragment), 1);
    Films = {'film','Film'};
    isFilm = ismember(Type, Films);
    Type(isFilm) = repmat({'film'}, sum(isFilm), 1);
    isOther = ~ismember(Type, [{'total'}, Fibres, Fragments, Films]);
    Type(isOther) = repmat({'other'}, sum(isOther), 1);
    abundance.Type = Type;
    plastics.(source).abundance = abundance;
end

% Omit variables not shared by all data tables -- these occur after the 
% Value column, and could be included later...
for i = 1:nsources
    source = sources{i};
    d = plastics.(source).abundance;
    d = d(:,1:find(strcmp(d.Properties.VariableNames, 'Value')));
    plastics.(source).abundance = d;
end

% Combine all plastic abundance data into a single table
abundance = structfun(@(z) z.abundance, plastics, 'UniformOutput', false);
abundance = struct2cell(abundance);
abundance = vertcat(abundance{:});

sources = unique(abundance.Source, 'stable');
nsources = length(sources);

% Save full table (and separate tables for different measurements
% variables)
saveData = true; % Save data table?
switch saveData, case true
    filename = 'plastic_quantity.csv';
    writetable(abundance, fullfile(baseDirectory, 'data/plastic_quantity', filename))
end
% switch Variables
%     case 'concentration', filename = 'plastic_conc.csv';
%     case 'density', filename = 'plastic_density.csv';
%     case 'mass density', filename = 'plastic_mass_density.csv';
%     case 'all', filename = 'plastic_quantity.csv';
% end



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
        pu = useUnit{i};
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
filename = 'plastic_mapPlot_SouthernOcean_Waller2017.png';
filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
exportgraphics(plt, filepath)

