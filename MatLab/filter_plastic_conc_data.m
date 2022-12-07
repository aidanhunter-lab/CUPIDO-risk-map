function [out, sources, nsources] = filter_plastic_conc_data(dat, varargin)
% Filters plastic abundance data sets and combines into a single table.
% Input 'dat' is a struct of data sets produced from
% load_plastic_conc_data.m

%% Optional arguments
extractVarargin(varargin)

if ~exist('onScreenDisplay', 'var')
    onScreenDisplay = false;
end
if ~exist('saveData', 'var')
    saveData = false;
end
if ~exist('outputFilename', 'var')
    filename = 'plastic_quantity.csv';
else
    filename = outputFilename;
end

if ~exist('returnOnlyMicroplasticMeasures', 'var')
    returnOnlyMicroplasticMeasures = false;
end

%% Directories
project = 'CUPIDO-risk-map';
thisFile = which('filter_plastic_conc_data.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);

%% Data sources
sources = fieldnames(dat);
nsources = length(sources);

%% Measurement units used by each source?
for i = 1:nsources
    source = sources{i};
    d = dat.(source).abundance;
    units.(source) = unique(d.Unit, 'stable')';
end
switch onScreenDisplay, case true, disp(units); end

% Measurement units may depend on sample type (sediment samples may be
% different).
uSampleTypes = [];
for i = 1:nsources
    source = sources{i};
    d = dat.(source).abundance;
    uSampleTypes = unique([uSampleTypes; d.SampleType], 'stable');
    sampleTypes.(source) = unique(d.SampleType, 'stable')';
end

% Volumetric concentration 'pieces/m^3' is most common unit for seawater 
% samples; some data use area density 'pieces/km^2', some use mass density
% 'g/km^2'.
varOpts_water = {'concentration', 'massConcentration', 'density', 'massDensity'}; %, 'all'}; % filtering options
varOpts_sediment = {'concentration', 'density'};

Units.water.concentration = {'pieces/m3', 'items/m3', 'particles/m3', 'number/m3'};
Units.water.massConcentration = {'g/m3'};
Units.water.density = {'pieces/km2', 'items/km2', 'number/km2'};
Units.water.massDensity = {'g/km2'};

Units.sediment.concentration = {'pieces/g'};
Units.sediment.density = {'pieces/m2'};

Units_water = struct2cell(Units.water);
useUnit_water = cellfun(@(z) z(1), Units_water);
Units_sediment = struct2cell(Units.sediment);
useUnit_sediment = cellfun(@(z) z(1), Units_sediment);

% Include a Variable column in each data set to indicate measurement type,
% omitting rows with units not listed in Units
for i = 1:nsources
    source = sources{i};
    sampleType = sampleTypes.(source);
    if contains(sampleType, 'water')
        Units_ = Units_water;
        useUnit = useUnit_water;
        varOpts = varOpts_water;
    end
    if contains(sampleType, 'sediment')
        Units_ = Units_sediment;
        useUnit = useUnit_sediment;
        varOpts = varOpts_sediment;
    end
    d = dat.(source).abundance;
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
    dat.(source).abundance = d;
end

%% Regularise the data and combine into a single table
for i = 1:nsources
    source = sources{i};
%     sampleType = sampleTypes.(source);
    d = dat.(source).abundance;
    fields = d.Properties.VariableNames;
    % Regularise the Station variable
    j = ismember(fields, 'station');
    if any(j), fields{j} = 'Station'; end
    j = ismember('Station', fields);
    if ~j
        k = strcmp(fields, 'Sample') | strcmp(fields, 'sample') | ... 
            strcmp(fields, 'SampleID') | strcmp(fields, 'MUC_ID');
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
    j = ismember(fields, {'Depth', 'depth', 'Dep', 'dep', 'Depth_m', 'depth_m'});
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
    d = movevars(d, 'SampleType', 'Before', 'Station');
    dat.(source).abundance = d;
end

% Only use microplastic measures?
switch returnOnlyMicroplasticMeasures, case true
    for i = 1:nsources
        source = sources{i};
        d = dat.(source).abundance;
        fields = d.Properties.VariableNames;
        j = ismember(fields, {'Size','size'});
        if ~any(j), continue; end
        j = ismember(d{:,j}, {'Micro', 'micro'});
        dat.(source).abundance = d(j,:);
    end
end

% Include a "Source" column in each abundance data set
for i = 1:nsources
    source = sources{i};
    d = dat.(source).abundance;
    source_ = source;
    j = strfind(source, '_');
    if ~isempty(j)
        source_(j) = ' ';
        source_ = [source_(1:end-4) '(' source_(end-3:end) ')'];
    end
    d.Source = repmat({source_}, height(d), 1);
    d = movevars(d, 'Source', 'Before', 'Station');
    dat.(source).abundance = d;
end

% Some data sources do not report the measurement means, but only the min 
% and max measurements at each station. This is usually because the data
% were reported as within some range [min, max]. Let's estimate mean values
% for those data where only min/max is reported simply by taking midpoints.

% 'M NOT KEEN ON THIS... SHOULD PROBABLY FIND A BETTER WAY TO DISPLAY THE
% DATA AS IT IS REPORTED BY EACH SOURCE

% useStat = 'max';
useStat = 'mean';
for i = 1:nsources
    source = sources{i};
    d = dat.(source).abundance;
    stats = unique(d.Stat, 'stable');
    useStatExists = ismember(useStat, stats);
    if useStatExists
        d = d(strcmp(d.Stat, useStat),:);
        dat.(source).abundance = d;
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
        dat.(source).abundance = d;
    end
end % The abundance data should now contain a single measurement per station

% Include info on plastic type (particle, fibre, etc).
for i = 1:nsources
    % Combine info on plastic category into the abundance data
    source = sources{i};
    abundance = dat.(source).abundance;
    categories = dat.(source).categories;
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
%             x_ = x.(Type)(:,{'Station', 'Unit', 'Value'});
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
    dat.(source).abundance = vertcat(abundance{:});
end

% Regularise the naming convention for plastic types
for i = 1:nsources
    source = sources{i};
    abundance = dat.(source).abundance;
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
    dat.(source).abundance = abundance;
end

% Some sample dates are not specific, but rather give an interval using 2
% columns. Therefore all Date variables need 2 columns
for i = 1:nsources
    source = sources{i};
    d = dat.(source).abundance;
    if size(d.Date, 2) == 1
        d.Date = [d.Date, NaT(height(d), 1)];
        dat.(source).abundance.Date = d.Date;
    end
end

% Omit variables not shared by all data tables -- these occur after the 
% Value column, and could be included later...
for i = 1:nsources
    source = sources{i};
    d = dat.(source).abundance;
    d = d(:,1:find(strcmp(d.Properties.VariableNames, 'Value')));
    dat.(source).abundance = d;
end

% Combine all plastic abundance data into a single table
out = structfun(@(z) z.abundance, dat, 'UniformOutput', false);
out = struct2cell(out);
out = vertcat(out{:});

sources = unique(out.Source, 'stable');
nsources = length(sources);

% Save full table (and separate tables for different measurements
% variables)
switch saveData, case true
    writetable(out, fullfile(baseDirectory, 'data/plastic_quantity', filename))
end

