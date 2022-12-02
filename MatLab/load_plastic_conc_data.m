function Data = load_plastic_conc_data(varargin)
% Load multiple sources of data on plastic concentration from Southern
% Ocean and store in a single struct

extractVarargin(varargin)

%% Set directories
if ~exist('project', 'var')
    project = 'CUPIDO-risk-map';
end
thisFile = which('load_plastic_conc_data');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
if ~exist('relDataDir', 'var')
    relDataDir = 'data/plastic_quantity';
end
dataDirectory = fullfile(baseDirectory, relDataDir);
addpath(genpath(dataDirectory));
addpath(genpath(fullfile(baseDirectory, 'MatLab')))


%% Relevant studies
%% Cincinelli 2017
source = 'Cincinelli_2017';
filepath = fullfile(dataDirectory, source);
filename = 'sample station locations_table 1.csv';
stations = readtable(fullfile(filepath, filename));

filename = 'plastic abundance at stations_figure 1.csv';
abundance = readtable(fullfile(filepath, filename));

filename = 'plastic categories.csv';
categories = readtable(fullfile(filepath, filename));

filename = 'plastic polymers.csv';
polymers = readtable(fullfile(filepath, filename));

% Total plastics
dat_tot = join(abundance, stations);
% dat_tot.Type = repmat({'total'}, height(dat_tot), 1);

dat_tot.Date = repmat(datetime({'01-Jan-2010', '28-Feb-2010'}), height(dat_tot), 1);
dat_tot.SampleType = repmat({'seawater'}, height(dat_tot), 1);
dat_tot = movevars(dat_tot, 'SampleType', 'Before', 'Sample');

% % Plastic categories (fibre, fragment, other)
% dat_cat = dat_tot; dat_cat(1:end,:) = [];
% Types = unique(categories.Type, 'stable');
% nTypes = length(Types);
% for i = min(dat_tot.Sample):max(dat_tot.Sample)
%     tot = dat_tot(dat_tot.Sample == i,:);
%     for j = 1:nTypes
%         p = 0.01 * categories.Value(strcmp(categories.Type, Types{j}) & strcmp(categories.Measure, 'Mean'));
%         d = tot;
%         d.Value = p .* d.Value;
%         d.Type = repmat(Types(j), height(d), 1);
%         dat_cat = vertcat(dat_cat, d);
%     end
% end
% 
% % Polymers
% dat_pol = dat_tot; dat_pol(1:end,:) = [];
% Polymers = unique(polymers.Polymer, 'stable');
% nPolymers = length(Polymers);
% for i = min(dat_tot.Sample):max(dat_tot.Sample)
%     tot = dat_tot(dat_tot.Sample == i,:);
%     for j = 1:nPolymers
%         p = 0.01 * polymers.RelativeAbundance(strcmp(polymers.Polymer, Polymers{j}));
%         d = tot;
%         d.Value = p .* d.Value;
%         d.Type = repmat(Polymers(j), height(d), 1);
%         dat_pol = vertcat(dat_pol, d);
%     end
% end

polymers.Properties.VariableNames{contains(polymers.Properties.VariableNames, 'Abundance')} = 'Value';

% Store all output in a struct
Cincinelli_2017.abundance = dat_tot;
Cincinelli_2017.polymers = polymers;
Cincinelli_2017.categories = categories;

Data.Cincinelli_2017 = Cincinelli_2017;

%% Isobe 2017

source = 'Isobe_2017';
filepath = fullfile(dataDirectory, source);
filename = 'plastic abundance_table 2.csv';
abundance = readtable(fullfile(filepath, filename));

% Lat-lon coords not provided in Isobe 2017, so fill in from inspection of
% figure 1.
coords = [107.5, -64; 117, -62.5; 128, -57.5; 134, -54.5; 142, -47.75];

abundance.Longitude = nan(height(abundance), 1);
abundance.Latitude = nan(height(abundance), 1);
stations = unique(abundance.Station, 'stable');
for i = 1:length(stations)
    abundance.Longitude(abundance.Station == stations(i)) = coords(i,1);
    abundance.Latitude(abundance.Station == stations(i)) = coords(i,2);
end

% All measurements are fragments
categories = removevars(abundance, {'Depth', 'Date', 'Variable', 'Value', 'Unit', 'Longitude', 'Latitude'});
categories = unique(categories);
categories.Value = repmat(100, height(categories), 1);
categories.Unit = repmat({'percent'}, height(categories), 1);
abundance.Type = [];
abundance.SampleType = repmat({'seawater'}, height(abundance), 1);
abundance = movevars(abundance, 'SampleType', 'Before', 'Station');

Isobe_2017.abundance = abundance;
Isobe_2017.categories = categories;

Data.Isobe_2017 = Isobe_2017;

%% Lacerda 2019

source = 'Lacerda_2019';
filepath = fullfile(dataDirectory, source);

filename = 'station locations.csv';
stations = readtable(fullfile(filepath, filename));

filename = 'plastic concentration_table S1.csv';
conc = readtable(fullfile(filepath, filename));

filename = 'plastic abundance by size class.csv';
abn_size = readtable(fullfile(filepath, filename));

filename = 'plastic abundance by polymer.csv';
abn_poly = readtable(fullfile(filepath, filename));

filename = 'plastic abundance by colour.csv';
abn_col = readtable(fullfile(filepath, filename));

filename = 'plastic abundance by type and size range.csv';
abn_type = readtable(fullfile(filepath, filename));

% Total plastics
conc.Properties.VariableNames{strcmp(conc.Properties.VariableNames, 'Sample')} = 'Station';
dat_tot = join(conc, stations);

abn_size.Properties.VariableNames{strcmp(abn_size.Properties.VariableNames, 'Abundance')} = 'Value';
abn_poly.Properties.VariableNames{strcmp(abn_poly.Properties.VariableNames, 'Abundance')} = 'Value';
abn_col.Properties.VariableNames{strcmp(abn_col.Properties.VariableNames, 'Abundance')} = 'Value';
abn_type.Properties.VariableNames{contains(abn_type.Properties.VariableNames, 'Abundance')} = 'Value';
abn_type.Unit = repmat({'percent'}, height(abn_type), 1);
abn_type = movevars(abn_type, 'Value', 'Before', 'Unit');
abn_type = movevars(abn_type, 'SizeRange', 'Before', 'Value');

dat_tot.Date = repmat(datetime({'01-Feb-2017', '28-Feb-2017'}), height(dat_tot), 1);
dat_tot.SampleType = repmat({'seawater'}, height(dat_tot), 1);
dat_tot = movevars(dat_tot, 'SampleType', 'Before', 'Station');

% Store all output in a struct
Lacerda_2019.abundance = dat_tot;
Lacerda_2019.size = abn_size;
Lacerda_2019.polymers = abn_poly;
Lacerda_2019.categories = abn_type;
Lacerda_2019.colour = abn_col;

Data.Lacerda_2019 = Lacerda_2019;

%% Jones-Williams 2020

source = 'Jones-Williams_2020';
filepath = fullfile(dataDirectory, source);

filename = 'sample metadata_table A1.csv';
meta = readtable(fullfile(filepath, filename));

filename = 'zooplankton concentration_table A2.csv';
opts = detectImportOptions(fullfile(filepath, filename));
opts = setvartype(opts, 'Station', 'char');
zoo = readtable(fullfile(filepath, filename), opts);

filename = 'plastic counts with size and type_table A3.csv';
counts = readtable(fullfile(filepath, filename));

filename = 'microplastic and microfibre concentration by station_table A6.csv';
conc = readtable(fullfile(filepath, filename));

stations = meta(strcmp(meta.SampleType, 'Water'),:);
stations.Properties.VariableNames{strcmp(stations.Properties.VariableNames, 'Sample')} = 'Station';
stations.Station = cell2mat(cellfun(@(z) str2double(z(9:end)), stations.Station, 'UniformOutput', false));

% To merge station info into conc table we need a single row per station,
% so average over sample start/end times
st = unique(stations.Station, 'stable');
nst = length(st);
d = stations; d(2:2:end,:) = [];
nc = size(d,2);
for i = 1:nst
    d_ = stations(stations.Station == st(i),:);
    for j = 1:nc
        if isnumeric(d_{:,j})
            d{i,j} = mean(d_{:,j});
        end
    end
end

dat = join(conc, d(:,{'Station','Date','Depth','Longitude','Latitude'}));

% As most other data sources here store info on plastic category in a
% separate table, it is useful to simply the abundance table by combining
% microplastics and microfibres and storing the proportions separately.
mplastic = dat(strcmp(dat.Type, 'microplastic'),:);
mfibre = dat(strcmp(dat.Type, 'microfibre'),:);

dat = mplastic;
dat.Value = mplastic.Value + mfibre.Value;
dat.Type = [];

dat.SampleType = repmat({'seawater'}, height(dat), 1);
dat = movevars(dat, 'SampleType', 'Before', 'Station');

categories = dat;
categories(3:3:height(categories),:) = [];
categories.Value = nan(height(categories), 1);
categories.Unit = repmat({'percent'}, height(categories), 1);
categories.Properties.VariableNames{'Variable'} = 'Type';
categories.Type = repmat({'microplastic'; 'microfibre'}, nst, 1);

for i = 1:nst
    mp = mplastic.Value(mplastic.Station == i & strcmp(mplastic.Variable, 'Concentration'),:);
    mf = mfibre.Value(mfibre.Station == i & strcmp(mfibre.Variable, 'Concentration'),:);
    categories.Value(categories.Station == i & strcmp(categories.Type, 'microplastic')) = 100 * mp / (mp + mf);
    categories.Value(categories.Station == i & strcmp(categories.Type, 'microfibre')) = 100 * mf / (mp + mf);
end

Jones_Williams_2020.metadata = meta;
Jones_Williams_2020.rawdata = counts;
Jones_Williams_2020.abundance = dat;
Jones_Williams_2020.categories = categories;
Jones_Williams_2020.zooplankton = zoo;

Data.Jones_Williams_2020 = Jones_Williams_2020;


%% Zhang 2022
% Really good, useful paper, but the data is not tabled so I'll need to
% extract it from their written results and figures... This  will take a
% while but it's worth it for samples from the eastern side of Antarctica

source = 'Zhang_2022';
filepath = fullfile(dataDirectory, source);

filename = 'plastic abundance at station_figure 2.csv';
stations = readtable(fullfile(filepath, filename));

filename = 'plastic abundance stats.csv';
stats = readtable(fullfile(filepath, filename));

filename = 'plastic colour.csv';
colours = readtable(fullfile(filepath, filename));

filename = 'plastic sizes.csv';
sizes = readtable(fullfile(filepath, filename));

filename = 'plastic type.csv';
types = readtable(fullfile(filepath, filename));

filename = 'polymers.csv';
polymers = readtable(fullfile(filepath, filename));

types.Properties.VariableNames(:,1:2) = {'Depth','Type'};
types.Type(strcmp(types.Type, 'Line/Fibre')) = {'Fibre'};

stations.Date = repmat(datetime({'01-Dec-2017', '31-Jan-2018'}), height(stations), 1);
stations.SampleType = repmat({'seawater'}, height(stations), 1);
stations = movevars(stations, 'SampleType', 'Before', 'Depth');

Zhang_2022.abundance = stations;
Zhang_2022.abundance_stats = stats;
Zhang_2022.size = sizes;
Zhang_2022.polymers = polymers;
Zhang_2022.size = sizes;
Zhang_2022.categories = types;
Zhang_2022.colour = colours;

Data.Zhang_2022 = Zhang_2022;


%% Buckingham 2022

source = 'Buckingham_2022';
filepath = fullfile(dataDirectory, source);
filename = 'SG_data_JBuckingham_v2.csv';
abundance = readtable(fullfile(filepath, filename));

abundance.Properties.VariableNames{'Measurement'} = 'Variable';
abundance.Properties.VariableNames{'Sample_Type'} = 'SampleType';
abundance.Properties.VariableNames{'Sample_Method'} = 'Sample_Method';
abundance.Properties.VariableNames{'Size_Class'} = 'SizeClass';

Site = unique(abundance.Site, 'stable');
Station = (1:length(Site))';
Stations = table(Station, Site);

size_dat = abundance(~strcmp(abundance.SizeClass, 'all'),:);
abundance(strcmp(abundance.SizeClass, 'all'),:);

abundance = join(abundance, Stations);
abundance = movevars(abundance, 'Station', 'Before', 'Site');

size_dat = join(size_dat, Stations);
size_dat = movevars(size_dat, 'Station', 'Before', 'Site');

% All measurements are fragments -- DOUBLE CHECK THIS
Type = {'particle'};
Value = 100;
Unit = {'percent'};
categories = table(Type, Value, Unit);

Buckingham_2022.abundance = abundance;
Buckingham_2022.size = size_dat;
Buckingham_2022.categories = categories;

Data.Buckingham_2022 = Buckingham_2022;


%% Waller 2017 providews a list of plastic samples with lat-lon coords.
source = 'Waller_2017';
filepath = fullfile(dataDirectory, source);
filename = 'tableS1.csv';
plastics = readtable(fullfile(filepath, filename));
plastics = renamevars(plastics, {'Type', 'Long', 'Lat', 'Quantity', 'Measure', 'Location_1'}, ...
    {'Size', 'Longitude', 'Latitude', 'Value' ,'Unit', 'Class'});

sources = unique(plastics.Source, 'stable');

%% Eriksen 2014

ind = strcmp(plastics.Source, sources{contains(sources, 'Eriksen')});
dat = plastics(ind,:);

dat = removevars(dat, {'Class'});
dat.Station = (1:height(dat))';

dat = movevars(dat, {'Station', 'Size'}, 'Before', 'Location');
dat = movevars(dat, {'Value', 'Unit'}, 'After', 'Location');
dat = movevars(dat, {'Longitude', 'Latitude'}, 'After', 'Unit');

dat.Unit(strcmp(dat.Unit, 'g km-2')) = {'g/km2'};
dat.Unit(strcmp(dat.Unit, 'Particles km-2')) = {'pieces/km2'};
dat.Size = strrep(dat.Size, 'M', 'm');

dat.SampleType = repmat({'seawater'}, height(dat), 1);
dat = movevars(dat, 'SampleType', 'Before', 'Station');

Eriksen_2014.abundance = dat;

% These are measurements of particles
Type = {'particle'};
Value = 100;
Unit = {'percent'};
Eriksen_2014.categories = table(Type, Value, Unit);

Data.Eriksen_2014 = Eriksen_2014;

%% Cozar 2014

ind = strcmp(plastics.Source, sources{contains(sources, 'Cózar')});
dat = plastics(ind,:);

dat = removevars(dat, {'Class'});
dat.Station = (1:height(dat))';

dat = movevars(dat, {'Station', 'Size'}, 'Before', 'Location');
dat = movevars(dat, {'Value', 'Unit'}, 'After', 'Location');
dat = movevars(dat, {'Longitude', 'Latitude'}, 'After', 'Unit');

dat.Unit(strcmp(dat.Unit, 'g km-2')) = {'g/km2'};
dat.Size = strrep(dat.Size, 'M', 'm');

dat.Date = repmat(datetime({'01-Dec-2012', '28-Feb-2013'}), height(dat), 1);
dat.SampleType = repmat({'seawater'}, height(dat), 1);
dat = movevars(dat, 'SampleType', 'Before', 'Station');

% Are these measurements of particles? check reference paper
Type = {'particle'};
Value = 100;
Unit = {'percent'};

Cozar_2014.abundance = dat;
Cozar_2014.categories = table(Type, Value, Unit);

Data.Cozar_2014 = Cozar_2014;



%% adventurescience.org

% I can extract adventurescience.org measures from data table provided by
% Waller, or I can load the values I extracted myself. Let's use my more
% recent extraction, although this lacks some of information provided in
% Waller's table.

source = 'AdventureScience';
filepath = fullfile(dataDirectory, source);

filename = 'plastic samples.csv';
abundance = readtable(fullfile(filepath, filename));

abundance.Sampling_Method = repmat({'Sample bottle'}, height(abundance), 1);
abundance.Depth = repmat({'<1m'}, height(abundance), 1);

abundance.SampleType(strcmp(abundance.Type, 'marine')) = {'seawater'};
abundance.SampleType(strcmp(abundance.Type, 'freshwater')) = {'freshwater'};
abundance.Type = [];
abundance = movevars(abundance, 'SampleType', 'Before', 'Sample');

Type = {'particle'};
Value = 100;
Unit = {'percent'};
categories = table(Type, Value, Unit);

adventurescience.abundance = abundance;
adventurescience.categories = categories;

Data.adventurescience = adventurescience;


% See commented code below for extracting adventurescience data from Waller's table
% ind = strcmp(plastics.Source, sources{contains(sources, 'adventure')});
% dat = plastics(ind,:);
% 
% dat = removevars(dat, {'Class'});
% dat.Station = (1:height(dat))';
% 
% dat = movevars(dat, {'Station', 'Size'}, 'Before', 'Location');
% dat = movevars(dat, {'Value', 'Unit'}, 'After', 'Location');
% dat = movevars(dat, {'Longitude', 'Latitude'}, 'After', 'Unit');
% 
% dat.Value = 1000 .* dat.Value; % convert particles/L -> particles/m^3
% dat.Unit = repmat({'particles/m3'}, height(dat), 1);
% dat.Size = strrep(dat.Size, 'M', 'm');
% 
% adventurescience.abundance = dat;
% 
% Type = {'particle'};
% Value = 100;
% Unit = {'percent'};
% adventurescience.categories = table(Type, Value, Unit);
% 
% Data.adventurescience = adventurescience;


%% Cunningham 2020 (sediment)

source = 'Cunningham_2020';
filepath = fullfile(dataDirectory, source);

filename = 'microplastic abundance_tableS1.csv';
abundance = readtable(fullfile(filepath, filename));

filename = 'microplastic properties_tableS2.csv';
properties = readtable(fullfile(filepath, filename));

filename = 'relative abundance of polymers_figure 2C.csv';
rel_abundance = readtable(fullfile(filepath, filename));

filename = 'sample locations_table1.csv';
stations = readtable(fullfile(filepath, filename));

filename = 'plastic type.csv';
categories = readtable(fullfile(filepath, filename));


abundance.Properties.VariableNames{'Measure'} = 'Replicate';
abundance.Unit = repmat({'pieces/g'}, height(abundance), 1);

stations.Year = nan(height(stations), 1);
stations.Year(contains(stations.Core, 'AP')) = 2017;
stations.Year(contains(stations.Core, {'SS', 'SG'}) & ...
    ~ismember(stations.MUC_ID, {'1-1', '7-2'})) = 2019;
stations.Year(ismember(stations.MUC_ID, {'1-1', '7-2'})) = 2019;

stations.Depth_m = arrayfun(@(z) {num2str(z)}, stations.Depth_m);

abundance = join(abundance, stations); % merge lat-lon and depth info into main table
abundance = movevars(abundance, {'Year', 'Longitude', 'Latitude', 'Depth_m'}, 'After', 'MUC_ID');
abundance = movevars(abundance, 'Unit', 'After', 'Mean');
abundance.SampleType = repmat({'sediment'}, height(abundance), 1);
abundance = movevars(abundance, 'SampleType', 'Before', 'Core');

rel_abundance.Properties.VariableNames([1,3]) = {'Plastic', 'Value'};

Cunningham_2020.abundance = abundance;
Cunningham_2020.polymers = rel_abundance;
Cunningham_2020.categories = categories;
Cunningham_2020.properties = properties;

% I NEED TO GET THE SAMPLE DATES FOR THESE SEDIMENT DATA -- IT'S NOT
% ENTIRELY CLEAR FROM THE PAPER SO I SHOULD EMAIL CUNNINGHAM TO ASK.
% HOWEVER, IT DOES LOOK AS THOUGH THE ANTARCTIC PENINSULA SAMPLES ARE FROM
% 2017 (JR17003A), THE SOUTH SANDWICH ISLANDS AND SOUTH GEORGIA SAMPLES ARE
% FROM 2019 (PS119 AND M134). NOTE THE DIFFERENCE IN SAMPLE ID AT THE
% BOTTOM OF TABLE 1 -- THIS MAY INDICATE THE DIFFERENT CRUISES IN 2019.

Data.Cunningham_2020 = Cunningham_2020;


%% Munari 2017 (sediment)

source = 'Munari_2017';
filepath = fullfile(dataDirectory, source);

filename = 'plastic abundance fractioned by size_table 4.csv';
abundance_size = readtable(fullfile(filepath, filename));

filename = 'plastic abundance fractioned by type and location_table 2.csv';
abundance_type = readtable(fullfile(filepath, filename));

filename = 'sample coords and depth_table 1.csv';
stations = readtable(fullfile(filepath, filename));

% convert coords from degree-minute into decimal degrees
lon = stations.Longtide_E;
lat = stations.Latitude_S;
lon = cell2mat(cellfun(@(z) str2double(z(1:3)) + str2double(z(5:end)) / 60, ...
    lon, 'UniformOutput', false));
lat = cell2mat(cellfun(@(z) -(str2double(z(1:3)) + str2double(z(5:end)) / 60), ...
    lat, 'UniformOutput', false));
stations.Longitude = lon;
stations.Latitude = lat;
stations = removevars(stations, {'Longtide_E','Latitude_S'});

sampleTime = datetime(2015,1,15); % samples collected in Jan 2015
stations.Date = repmat(sampleTime, height(stations), 1);

abundance_type.Properties.VariableNames{'Site'} = 'Station';
abundance_size.Properties.VariableNames{'Site'} = 'Station';

abundance_type = join(abundance_type, stations);
abundance_size = join(abundance_size, stations);

abundance_type = movevars(abundance_type, {'Longitude','Latitude','Depth_m', 'Date'}, 'After', 'Station');
abundance_size = movevars(abundance_size, {'Longitude','Latitude','Depth_m', 'Date'}, 'After', 'Station');

abundance = abundance_type;

m = reshape(abundance.Value, 3, []);
m = sum(m);
abundance = removevars(abundance, {'Plastic', 'Value'});
abundance = unique(abundance, 'stable');
abundance.Value = m(:);
abundance = movevars(abundance, 'Value', 'Before', 'Unit');

abundance.SampleType = repmat({'sediment'}, height(abundance), 1);
abundance = movevars(abundance, 'SampleType', 'Before', 'Location');

abundance_type.Properties.VariableNames{'Plastic'} = 'Type';

Station = unique(abundance.Station, 'stable');
snew = (1:length(Station))';
x = table(Station, snew);

abundance = join(abundance, x);
abundance = removevars(abundance, 'Station');
abundance.Properties.VariableNames{'snew'} = 'Station';
abundance = movevars(abundance, 'Station', 'After', 'Location');

abundance_type = join(abundance_type, x);
abundance_type = removevars(abundance_type, 'Station');
abundance_type.Properties.VariableNames{'snew'} = 'Station';
abundance_type = movevars(abundance_type, 'Station', 'After', 'Location');

abundance_size = join(abundance_size, x);
abundance_size = removevars(abundance_size, 'Station');
abundance_size.Properties.VariableNames{'snew'} = 'Station';
abundance_size = movevars(abundance_size, 'Station', 'After', 'Location');

abundance.Depth_m = arrayfun(@(z) {num2str(z)}, abundance.Depth_m);

x.Station = [];
x.Station = x.snew;
x.snew = [];

Types = unique(abundance_type.Type);
nt = length(Types);
Station = reshape(repmat(x.Station', nt, 1), [], 1);
AT = table(Station);
for i = 1:height(x)
    for j = 1:nt
        k = strcmp(abundance_type.Stat, 'Mean') & abundance_type.Station == i & strcmp(abundance_type.Type, Types{j});
        v = abundance_type.Value(k);
        k0 = strcmp(abundance.Stat, 'Mean') & abundance.Station == i;
        v0 = abundance.Value(k0);
        p = 100 * v / v0;
        k = find(AT.Station == i);
        k = k(j);
        AT.Type(k) = Types(j);
        AT.Value(k) = p;
    end
end
AT.Unit = repmat({'percent'}, height(AT), 1);

Sizes = unique(abundance_size.Size);
nt = length(Sizes);
Station = reshape(repmat(x.Station', nt, 1), [], 1);
ST = table(Station);
for i = 1:height(x)
    for j = 1:nt
        k = strcmp(abundance_size.Variable, 'Mean') & abundance_size.Station == i & strcmp(abundance_size.Size, Sizes{j});
        v = abundance_size.Value(k);
        k0 = strcmp(abundance.Stat, 'Mean') & abundance.Station == i;
        v0 = abundance.Value(k0);
        p = 100 * v / v0;
        k = find(ST.Station == i);
        k = k(j);
        ST.Size(k) = Sizes(j);
        ST.Value(k) = p;
    end
end
ST.Unit = repmat({'percent'}, height(ST), 1);


Munari_2017.abundance = abundance;
Munari_2017.size = ST;
Munari_2017.categories = AT;

Data.Munari_2017 = Munari_2017;















% Suaria 2020


% Studies cited in Waller 2017
% Eriksen 2014 (Gloabl study -- possibly poor data for Southern Ocean)
% Isobe 2017
% Cozar 2014 (Another global study -- looks a few samples from Antartic Peninsula)
% Barnes 2010

% Other notes...
% See Everaert (2018) for micropalstic risk assessment
% See Hardesty (2017) for micropalstic modelling review
% See Belcher (2019) for krill POC modelling (defo read this)
% See Bohdan for lit review on plastic abundance

