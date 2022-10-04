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

%% Load data

% Waller (2017) providews a list of plastic samples with lat-lon coords.
% Just use this for now -- I can include more later...
source = 'Waller_2017';
filepath = fullfile(dataDirectory, source);
filename = 'tableS1.csv';
plastics = readtable(fullfile(filepath, filename));

unique(plastics.Source, 'stable');

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
% Really good, useful paper, but the data is not tabled sp I'll need to
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

Zhang_2022.abundance = stations;
Zhang_2022.abundance_stats = stats;
Zhang_2022.size = sizes;
Zhang_2022.polymers = polymers;
Zhang_2022.size = sizes;
Zhang_2022.categories = types;
Zhang_2022.colour = colours;

Data.Zhang_2022 = Zhang_2022;



% Suaria 2020

% Cunningham 2020 (sediment)
% Munari 2017 (sediment)

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

% Isobe (2017) has method for converting surface MP concentration into
% total water concentration