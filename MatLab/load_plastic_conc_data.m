function [Data, DATA] = load_plastic_conc_data(varargin)
% Load multiple sources of data on plastic concentration from Southern
% Ocean and store all info in a single struct, Data, and abundance info in
% a table, DATA.

extractVarargin(varargin)

% Standardise the format of data from different sources so they may be
% combined into a single table.

% Data columns:
columnDescriptions.Source          = 'First column: identifies the study that data is sourced from';
columnTypes.Source                 = 'cellstr';
columnDescriptions.SampleType      = 'Sample medium -- seawater, sediment etc';
columnTypes.SampleType             = 'cellstr';
columnDescriptions.SampleGear      = 'Gear used for sampling -- bottles, nets, etc';
columnTypes.SampleGear             = 'cellstr';
columnDescriptions.LitterIDMethod  = 'Technique for identifying microplastic/litter -- FTIR, visual, etc';
columnTypes.LitterIDMethod         = 'cellstr';
columnDescriptions.LitterCategory  = 'Type of litter/particle -- plastic, metal, lithogenic, etc';
columnTypes.LitterCategory         = 'cellstr';
columnDescriptions.LitterScale     = 'Size class of litter -- nano, micro, meso, macro';
columnTypes.LitterScale            = 'cellstr';
columnDescriptions.PlasticForm     = 'Structure of plastic -- fragment, fibre, etc';
columnTypes.PlasticForm            = 'cellstr';
columnDescriptions.PlasticSize     = 'Measured size of plastic -- some numeric values, some categorical => store as character variable';
columnTypes.PlasticSize            = 'cellstr';
columnDescriptions.SampleAtStation = 'Logical variable: true if sample is from fixed station, false if sample derives from tow over long distance';
columnTypes.SampleAtStation        = 'logical';
columnDescriptions.Site            = 'Description of sample location -- character variable (empty values are allowed)';
columnTypes.Site                   = 'cellstr';
columnDescriptions.SiteCategory    = 'Categorical variable describing sample site, e.g., nearshore, offshore (empty values are allowed)';
columnTypes.SiteCategory           = 'cellstr';
columnDescriptions.SampleID        = 'Sample identifier -- unique for each station, event or tow';
columnTypes.SampleID               = 'cellstr';
columnDescriptions.Date            = 'Sample date: sometimes a precise date for each sample, otherwise a date range covering the sampling campaign';
columnTypes.Date                   = 'datetime'; % column width = 2
columnDescriptions.Longitude       = 'Coordinates (decimal degrees) -- take average values for long tows where SampleAtStation = false';
columnTypes.Longitude              = 'double';
columnDescriptions.Latitude        = 'Coordinates (decimal degrees) -- take average values for long tows where SampleAtStation = false';
columnTypes.Latitude               = 'double';
columnDescriptions.Longitude_start = 'Coordinates at start of tow -- applicable to SampleAtStation = false, NaN for SampleAtStation = true';
columnTypes.Longitude_start        = 'double';
columnDescriptions.Latitude_start  = 'Coordinates at start of tow -- applicable to SampleAtStation = false, NaN for SampleAtStation = true';
columnTypes.Latitude_start         = 'double';
columnDescriptions.Longitude_end   = 'Coordinates at start of tow -- applicable to SampleAtStation = false, NaN for SampleAtStation = true';
columnTypes.Longitude_end          = 'double';
columnDescriptions.Latitude_end    = 'Coordinates at start of tow -- applicable to SampleAtStation = false, NaN for SampleAtStation = true';
columnTypes.Latitude_end           = 'double';
columnDescriptions.Depth           = 'Sample depth -- character variable, a precise numeric value (in metres), otherwise a string variable';
columnTypes.Depth                  = 'cellstr';
columnDescriptions.Variable        = 'Measurement type -- concentration, density, etc';
columnTypes.Variable               = 'cellstr';
columnDescriptions.Statistic       = 'Measurement statistic -- categorical variable: most often raw but also mean, median, min/max';
columnTypes.Statistic              = 'cellstr';
columnDescriptions.Replicate       = 'Sample number within SampleID -- numeric value, usually 1 as single samples are most common';
columnTypes.Replicate              = 'double';
columnDescriptions.Value           = 'Numeric measurement value';
columnTypes.Value                  = 'double';
columnDescriptions.Unit            = 'Measurement unit corresponding to Value column -- use SI units';
columnTypes.Unit                   = 'cellstr';
columnDescriptions.Observation     = 'Additional comments -- store text descriptions of Value here, e.g., presence/absence data';
columnTypes.Observation            = 'cellstr';

allColumns = fieldnames(columnDescriptions);



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

polymers.Properties.VariableNames{contains(polymers.Properties.VariableNames, 'Abundance')} = 'Value';

% Merge data on plastic form
pc = unique(categories.Type);
dat_tot.PlasticForm = repmat({'all'}, height(dat_tot), 1);
dat_tot_ = dat_tot;
for i = 1:length(pc)
    d = dat_tot_;
    d.PlasticForm = repmat(pc(i), height(dat_tot_), 1);
    dv = categories.Value(strcmp(categories.Type, pc(i)) & strcmp(categories.Measure, 'Mean'));
    du = categories.Unit{strcmp(categories.Type, pc(i)) & strcmp(categories.Measure, 'Mean')};
    switch du
        case {'percent','Percent'}
            d.Value = 0.01 .* dv .* d.Value;
        case {'proportion','Proportion',''}
            d.Value = dv .* d.Value;
    end
    dat_tot = vertcat(dat_tot, d);
end


% Regularise abundance data
hd = height(dat_tot);
dat_tot.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
dat_tot.SampleGear = repmat({'vessel underway pump'}, hd, 1);
dat_tot.LitterIDMethod = repmat({'FTIR'}, hd, 1);
dat_tot.LitterCategory = repmat({'plastic'}, hd, 1);
dat_tot.LitterScale = repmat({'micro'}, hd, 1);
dat_tot.PlasticSize = cell(hd, 1); dat_tot.PlasticSize(:) = {''};
dat_tot.SampleAtStation = true(hd, 1);
dat_tot.Site = cell(hd, 1); dat_tot.Site(:) = {''};
dat_tot.SiteCategory = cell(hd, 1); dat_tot.SiteCategory(:) = {''};
dat_tot.SampleID = cellstr(num2str(dat_tot.Sample));
dat_tot.Longitude_start = nan(hd, 1);
dat_tot.Latitude_start = nan(hd, 1);
dat_tot.Longitude_end = nan(hd, 1);
dat_tot.Latitude_end = nan(hd, 1);
dat_tot.Variable = repmat({'concentration'}, hd, 1);
dat_tot.Statistic = dat_tot.Measure;
dat_tot.Replicate = ones(hd, 1);
dat_tot.Observation = cell(hd, 1); dat_tot.Observation(:) = {''};

dat_tot = dat_tot(:,allColumns);

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

% Size distribution (averaged over stations) of sampled plastics can be
% extracted from figure 3...


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

% Regularise abundance data
hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = repmat({'neutson net 0.35 mm'}, hd, 1);
abundance.LitterIDMethod = repmat({'FTIR'}, hd, 1);
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = repmat({'micro'}, hd, 1);
abundance.PlasticForm = repmat({'fragment'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = cell(hd, 1); abundance.Site(:) = {''};
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = cellstr(num2str(abundance.Station));
abundance.Date = [abundance.Date, NaT(hd, 1)];
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

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
% see figure 1 for station-specific distributions
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

% Merge data on plastic form
pc = unique(abn_type.Type);
dat_tot.PlasticForm = repmat({'all'}, height(dat_tot), 1);
dat_tot_ = dat_tot;
for i = 1:length(pc)
    d = dat_tot_;
    d.PlasticForm = repmat(pc(i), height(dat_tot_), 1);
    dv = abn_type.Value(strcmp(abn_type.Type, pc(i)));
    du = abn_type.Unit{strcmp(abn_type.Type, pc(i))};
    switch du
        case {'percent','Percent'}
            d.Value = 0.01 .* dv .* d.Value;
        case {'proportion','Proportion',''}
            d.Value = dv .* d.Value;
    end
    dat_tot = vertcat(dat_tot, d);
end

abundance = dat_tot;

hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = cell(hd, 1); abundance.SampleGear(:) = {''};
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = repmat({'micro and meso -- see size data'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = cell(hd, 1); abundance.Site(:) = {''};
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = cellstr(num2str(abundance.Station));
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);



% Store all output in a struct
Lacerda_2019.abundance = abundance;
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


% Merge data on plastic form
pc = unique(categories.Type);
dat.PlasticForm = repmat({'all'}, height(dat), 1);
dat_ = dat;
ns = unique(dat.Station);
for i = 1:length(pc)
    d = dat_;
    d.PlasticForm = repmat(pc(i), height(dat_), 1);
    for j = 1:length(ns)
        k = ns(j);
        dv = categories.Value(categories.Station == k & strcmp(categories.Type, pc{i}));
        du = categories.Unit{categories.Station == k & strcmp(categories.Type, pc{i})};
        switch du
            case {'percent','Percent'}
                d.Value(d.Station == k) = 0.01 .* dv .* d.Value(d.Station == k);
            case {'proportion','Proportion',''}
                d.Value(d.Station == k) = dv .* d.Value(d.Station == k);
        end
    end
    dat = vertcat(dat, d);
end

abundance = dat;

hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = cell(hd, 1); abundance.SampleGear(:) = {''};
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = repmat({'micro'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = cell(hd, 1); abundance.Site(:) = {''};
abundance.SiteCategory = abundance.Location;
abundance.SampleID = cellstr(num2str(abundance.Station));
abundance.Date = [abundance.Date, NaT(hd, 1)];
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

Jones_Williams_2020.metadata = meta;
Jones_Williams_2020.rawdata = counts;
Jones_Williams_2020.abundance = abundance;
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
% types.Type(strcmp(types.Type, 'Line/Fibre')) = {'Fibre'};

stations.Date = repmat(datetime({'01-Dec-2017', '31-Jan-2018'}), height(stations), 1);
stations.SampleType = repmat({'seawater'}, height(stations), 1);
stations = movevars(stations, 'SampleType', 'Before', 'Depth');

abundance = stations;

% Merge data on plastic form
pc = unique(types.Type);
abundance.PlasticForm = repmat({'all'}, height(abundance), 1);
abundance_ = abundance;
nd = unique(abundance.Depth);
for i = 1:length(pc)
    d = abundance_;
    d.PlasticForm = repmat(pc(i), height(abundance_), 1);
    for j = 1:length(nd)
        k = nd{j};
        dv = types.Value(strcmp(types.Depth, k) & strcmp(types.Type, pc{i}));
        du = types.Unit{strcmp(types.Depth, k) & strcmp(types.Type, pc{i})};
        switch du
            case {'percent','Percent'}
                d.Value(strcmp(d.Depth, k)) = 0.01 .* dv .* d.Value(strcmp(d.Depth, k));
            case {'proportion','Proportion',''}
                d.Value(strcmp(d.Depth, k)) = dv .* d.Value(strcmp(d.Depth, k));
        end
    end
    abundance = vertcat(abundance, d);
end


hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = cell(hd, 1); abundance.SampleGear(:) = {''};
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = repmat({'micro'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = abundance.Location;
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = abundance.SampleID;
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = abundance.Stat;
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

Zhang_2022.abundance = abundance;
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

micro = ismember(abundance.SizeClass, {'50-100 mum', '100-300 mum', '300-1000 mum'});
meso = ismember(abundance.SizeClass, {'1000-5000 mum'});
micromeso = ismember(abundance.SizeClass, {'all'});
abundance.LitterScale(micro) = {'micro'};
abundance.LitterScale(meso) = {'meso'};
abundance.LitterScale(micromeso) = {'micro and meso'};

size_dat = join(size_dat, Stations);
size_dat = movevars(size_dat, 'Station', 'Before', 'Site');

% All measurements are fragments
Type = {'particle'};
Value = 100;
Unit = {'percent'};
categories = table(Type, Value, Unit);


hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = abundance.Sample_Method;
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
% abundance.LitterScale = repmat({'micro'}, hd, 1);
abundance.PlasticForm = repmat({'particle'}, hd, 1);
abundance.PlasticSize = abundance.SizeClass;
abundance.SampleAtStation = true(hd, 1);
% abundance.Site = abundance.Location;
abundance.SiteCategory = abundance.Group;
abundance.SampleID = cellstr(num2str(abundance.Station));
abundance.Date = [abundance.Date, NaT(hd, 1)];
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

Buckingham_2022.abundance = abundance;
Buckingham_2022.size = size_dat;
Buckingham_2022.categories = categories;

Data.Buckingham_2022 = Buckingham_2022;


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

hd = height(abundance);
abundance.Source = repmat({source}, hd, 1);
abundance.SampleGear = abundance.Sampling_Method;
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = cell(hd, 1); abundance.LitterScale(:) = {''};
abundance.PlasticForm = repmat({'particle'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = cell(hd, 1); abundance.Site(:) = {''};
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = cellstr(num2str(abundance.Sample));
abundance.Date = [abundance.Date, NaT(hd, 1)];
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
% abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

AdventureScience.abundance = abundance;
AdventureScience.categories = categories;

Data.AdventureScience = AdventureScience;


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

dat.Variable(strcmp(dat.Unit, 'g/km2')) = {'massDensity'};
dat.Variable(strcmp(dat.Unit, 'pieces/km2')) = {'density'};

dat.Date = repmat(datetime({'01-Jan-2011', '31-Dec-2013'}), height(dat), 1);

abundance = dat;

hd = height(abundance);
abundance.Source = repmat({'Eriksen (2014)'}, hd, 1);
abundance.SampleGear = abundance.Sampling_Method;
abundance.LitterIDMethod = abundance.Identification_Method;
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = abundance.Size;
abundance.PlasticForm = repmat({'particle'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = abundance.Location;
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = cellstr(num2str(abundance.Station));
% abundance.Date = NaT(hd, 2); % IMPORTANT: FIND THE DATE FOR THESE SAMPLES
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
% abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

Eriksen_2014.abundance = abundance;

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

dat.Variable = repmat({'massDensity'}, height(dat), 1);

% Are these measurements of particles? check reference paper
Type = {'particle'};
Value = 100;
Unit = {'percent'};
categories = table(Type, Value, Unit);

abundance = dat;

hd = height(abundance);
abundance.Source = repmat({'Cozar (2014)'}, hd, 1);
abundance.SampleGear = abundance.Sampling_Method;
abundance.LitterIDMethod = abundance.Identification_Method;
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = abundance.Size;
abundance.PlasticForm = repmat({'plastic'}, hd, 1);
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = abundance.Location;
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = cellstr(num2str(abundance.Station));
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
% abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);


Cozar_2014.abundance = abundance;
Cozar_2014.categories = categories;

Data.Cozar_2014 = Cozar_2014;


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

% JR17003a cruise info 
% https://www.bodc.ac.uk/resources/inventories/cruise_inventory/reports/jr17003a.pdf
% PS119 cruise info
% https://www.tib.eu/en/suchen/id/awi:2f4972cb5868f267f79cd02562f3f9cbd05a15ef
% M134 cruise info
% https://www.bodc.ac.uk/resources/inventories/cruise_inventory/reports/meteor_m134.pdf

stations.Depth_m = arrayfun(@(z) {num2str(z)}, stations.Depth_m);

abundance.Properties.VariableNames{'Measure'} = 'Replicate';
abundance.Unit = repmat({'pieces/g'}, height(abundance), 1);

% stations.Year = nan(height(stations), 1);
% stations.Year(contains(stations.Core, 'AP')) = 2017;
% stations.Year(contains(stations.Core, {'SS', 'SG'}) & ...
%     ~ismember(stations.MUC_ID, {'1-1', '7-2'})) = 2019;
% stations.Year(ismember(stations.MUC_ID, {'1-1', '7-2'})) = 2019;

% stations.Depth_m = arrayfun(@(z) {num2str(z)}, stations.Depth_m);

abundance = join(abundance, stations); % merge lat-lon and depth info into main table
abundance = movevars(abundance, 'Cruise', 'Before', 'Core');
abundance = movevars(abundance, {'Date', 'Longitude', 'Latitude', 'Depth_m'}, 'After', 'MUC_ID');
abundance = movevars(abundance, 'Unit', 'After', 'Mean');
abundance.SampleType = repmat({'sediment'}, height(abundance), 1);
abundance = movevars(abundance, 'SampleType', 'Before', 'Cruise');

abundance.Variable = repmat({'concentration'}, height(abundance), 1);
abundance.Replicate = strrep(abundance.Replicate, 'sample ', '');
abundance.Replicate = cellfun(@(z) str2double(z), abundance.Replicate);

rel_abundance.Properties.VariableNames([1,3]) = {'Plastic', 'Value'};

% Merge data on plastic form
pc = unique(categories.Type);
abundance.PlasticForm = repmat({'all'}, height(abundance), 1);
abundance_ = abundance;
for i = 1:length(pc)
    d = abundance_;
    d.PlasticForm = repmat(pc(i), height(abundance_), 1);
    dv = categories.Value(strcmp(categories.Type, pc{i}));
    du = categories.Unit{strcmp(categories.Type, pc{i})};
    switch du
        case {'percent','Percent'}
            d.Value = 0.01 .* dv .* d.Value;
        case {'proportion','Proportion',''}
            d.Value = dv .* d.Value;
    end
    abundance = vertcat(abundance, d);
end

hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = repmat({'OKTOPUS multicore'}, hd, 1);
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = cell(hd, 1); abundance.LitterScale(:) = {''};
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = cell(hd, 1); abundance.Site(:) = {''};
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = abundance.MUC_ID; % could also use Core
abundance.Date = [abundance.Date, NaT(hd, 1)];
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Depth = abundance.Depth_m;
% abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = repmat({'raw'}, hd, 1);
% abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

Cunningham_2020.abundance = abundance;
Cunningham_2020.polymers = rel_abundance;
Cunningham_2020.categories = categories;
Cunningham_2020.properties = properties;

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


Mean = reshape(abundance.Value(strcmp(abundance.Stat, 'Mean')), 3, []);
StdDev = reshape(abundance.Value(strcmp(abundance.Stat, 'StdDev')), 3, []);
Mean = sum(Mean); Mean = Mean(:);
StdDev = sum(StdDev .^ 2) .^ 0.5; StdDev = StdDev(:);
abundance_ = removevars(abundance, {'Plastic', 'Value'});
abundance_ = unique(abundance_, 'stable');
abundance_.Value = [Mean; StdDev];
abundance_.Plastic = repmat({'all'}, height(abundance_), 1);
abundance = vertcat(abundance_, abundance);


abundance.SampleType = repmat({'sediment'}, height(abundance), 1);
abundance = movevars(abundance, 'SampleType', 'Before', 'Location');

abundance.Depth_m = cellstr(num2str(abundance.Depth_m));

abundance.Variable = repmat({'density'}, height(abundance), 1);


hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = repmat({'Van Veen grab'}, hd, 1);
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = repmat({'plastic'}, hd, 1);
abundance.LitterScale = repmat({'see size data'}, hd, 1);
abundance.PlasticForm = abundance.Plastic;
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
abundance.SampleAtStation = true(hd, 1);
abundance.Site = abundance.Location;
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = abundance.Station;
abundance.Date = [abundance.Date, NaT(hd, 1)];
abundance.Longitude_start = nan(hd, 1);
abundance.Latitude_start = nan(hd, 1);
abundance.Longitude_end = nan(hd, 1);
abundance.Latitude_end = nan(hd, 1);
abundance.Depth = abundance.Depth_m;
% abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = abundance.Stat;
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);


Munari_2017.abundance = abundance;
Munari_2017.size = abundance_size;
Munari_2017.categories = abundance_type;

Data.Munari_2017 = Munari_2017;


%% Suaria 2020


source = 'Suaria_2020';
filepath = fullfile(dataDirectory, source);

filename = 'All neuston samples table S1.csv';
abundance = readtable(fullfile(filepath, filename));

filename = 'microplastic categories.csv';
mp_cat = readtable(fullfile(filepath, filename));

filename = 'microplastic colour.csv';
mp_col = readtable(fullfile(filepath, filename));

filename = 'microplastic size.csv';
mp_size = readtable(fullfile(filepath, filename));

filename = 'polymer composition 0-25mm table 2.csv';
mp_com = readtable(fullfile(filepath, filename));

filename = 'Visual survey table S2.csv';
visual = readtable(fullfile(filepath, filename));

filename = 'Macrolitter metadata table S3.csv';
macro = readtable(fullfile(filepath, filename));


abundance.SampleType = repmat({'seawater'}, height(abundance), 1);
abundance.Depth = repmat({'surface'}, height(abundance), 1);

abundance.SampleMethod = repmat({'neuston net 200 mu m'}, height(abundance), 1);

% Adjust units from /L to /m3
j = contains(abundance.Unit, '/L');
abundance.Value(j) = 1e3 .* abundance.Value(j);
abundance.Unit(j) = strrep(abundance.Unit(j), '/L', '/m3');

% microplastic categories
mp_cat.Properties.VariableNames{'Category'} = 'Type';

% Convert longitude to degrees east (0,360)
w = visual.Longitude_start < 0;
visual.Longitude_start(w) = 180 + (180 + visual.Longitude_start(w));
w = visual.Longitude_end < 0;
visual.Longitude_end(w) = 180 + (180 + visual.Longitude_end(w));
w = macro.Longitude < 0;
macro.Longitude(w) = 180 + (180 + macro.Longitude(w));

% Merge sampleID info macro
macro.SampleID = cell(height(macro), 1); macro.SampleID(:) = {''};
v = visual(:,{'SampleID', 'Latitude_start', 'Longitude_start', 'Latitude_end', 'Longitude_end', 'Observation_point'});
v = unique(v, 'stable');
v.vdir = -(v.Latitude_start - v.Latitude_end) > 0; % indicate tow direction (true = north)
v.hdir = -(v.Longitude_start - v.Longitude_end) > 0; % indicate tow direction (true = east)
for i = 1:height(macro)
    m = macro(i,:);
    ilat = (v.vdir & v.Latitude_start <= m.Latitude & m.Latitude <= v.Latitude_end) | ...
        (~v.vdir & v.Latitude_end <= m.Latitude & m.Latitude <= v.Latitude_start);
    ilon = (v.hdir & v.Longitude_start <= m.Longitude & m.Longitude <= v.Longitude_end) | ...
        (~v.hdir & v.Longitude_end <= m.Longitude & m.Longitude <= v.Longitude_start);
%     io = strcmp(v.Observation_point, m.Observation_point);
    j = ilat & ilon; % & io;
    sj = sum(j);
    if sj == 0, match = 'no'; elseif sj == 1, match = 'yes'; elseif sj > 1, match = 'multiple'; end
    switch match
        case 'yes'
            macro.SampleID{i,:} = v.SampleID{j};
        case 'multiple'
            macro.SampleID{i,:} = 'multiple';
    end
end
macro = macro(~strcmp(macro.SampleID, 'multiple'),:);
macro = macro(~cellfun(@(z) isempty(z), macro.SampleID),:);

% Merge the microplastic and macro-litter abundance tables
visual.Depth = repmat({'surface'}, height(visual), 1);
visual.Station = visual.SampleID;
% visual.Latitude = nan(height(visual), 1);
% visual.Longitude = nan(height(visual), 1);
visual.Latitude = 0.5 .* (visual.Latitude_start + visual.Latitude_end); % tow midpoints position data on map
visual.Longitude = 0.5 .* (visual.Longitude_start + visual.Longitude_end);
visual.SampleType = repmat({'seawater'}, height(visual), 1);
visual.Measure = cell(height(visual), 1);
visual.Measure(strcmp('density', visual.Variable)) = {'Mean'};
visual.Measure(contains(visual.Variable, 'lower')) = {'Min'};
visual.Measure(contains(visual.Variable, 'upper')) = {'Max'};
visual.Variable(strcmp(visual.Variable, 'mass density_lower')) = {'massDensity'};
visual.Variable(strcmp(visual.Variable, 'mass density_upper')) = {'massDensity'};
abundance.Latitude_start = nan(height(abundance), 1);
abundance.Latitude_end = abundance.Latitude_start;
abundance.Longitude_start = abundance.Latitude_start;
abundance.Longitude_end = abundance.Latitude_start;
abundance.TowDistance_km = nan(height(abundance), 1);
abundance.Observation_point = cell(height(abundance), 1);
abundance.Group = repmat({'plastic'}, height(abundance), 1);
abundance.Measure = repmat({'Mean'}, height(abundance), 1);
abundance.Station = cellstr(num2str(abundance.Station));

abundance = abundance(:,[10 1:2 18 3 16:17 4:6 12:15 11 7 19 8:9]);
visual = visual(:,abundance.Properties.VariableNames);


% Merge data on plastic form -- microplastic
pc = unique(mp_cat.Type);
abundance.PlasticForm = repmat({'all'}, height(abundance), 1);
abundance_ = abundance;
for i = 1:length(pc)
    d = abundance_;
    d.PlasticForm = repmat(pc(i), height(abundance_), 1);
    dv = mp_cat.Value(strcmp(mp_cat.Type, pc{i}));
    du = mp_cat.Unit{strcmp(mp_cat.Type, pc{i})};
    switch du
        case {'percent','Percent'}
            d.Value = 0.01 .* dv .* d.Value;
        case {'proportion','Proportion',''}
            d.Value = dv .* d.Value;
    end
    abundance = vertcat(abundance, d);
end
% % Merge data on plastic form -- macrolitter
% pc = unique(macro.Category);
% abundance.PlasticForm = repmat({'all'}, height(abundance), 1);
% abundance_ = abundance;
% for i = 1:length(pc)
%     d = abundance_;
%     d.PlasticForm = repmat(pc(i), height(abundance_), 1);
%     dv = mp_cat.Value(strcmp(mp_cat.Type, pc{i}));
%     du = mp_cat.Unit{strcmp(mp_cat.Type, pc{i})};
%     switch du
%         case {'percent','Percent'}
%             d.Value = 0.01 .* dv .* d.Value;
%         case {'proportion','Proportion',''}
%             d.Value = dv .* d.Value;
%     end
%     abundance = vertcat(abundance, d);
% end

visual.PlasticForm = cell(height(visual), 1); visual.PlasticForm(:) = {''};

abundance = [abundance; visual]; % merge

abundance.Class(contains(abundance.Class, 'macro')) = {'macro'};
abundance.Class(contains(abundance.Class, 'micro')) = {'micro'};

abundance.SampleAtStation(isnan(abundance.TowDistance_km)) = true;
abundance.SampleAtStation(~isnan(abundance.TowDistance_km)) = false;

hd = height(abundance);
abundance.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
abundance.SampleGear = abundance.SampleMethod;
abundance.LitterIDMethod = cell(hd, 1); abundance.LitterIDMethod(:) = {''};
abundance.LitterCategory = abundance.Group;
abundance.LitterScale = abundance.Class;
abundance.PlasticSize = cell(hd, 1); abundance.PlasticSize(:) = {''};
% abundance.SampleAtStation = true(hd, 1);
abundance.Site = cell(hd, 1); abundance.Site(:) = {''};
abundance.SiteCategory = cell(hd, 1); abundance.SiteCategory(:) = {''};
abundance.SampleID = abundance.Station;
abundance.Date = [abundance.Date, NaT(hd, 1)];
% abundance.Variable = repmat({'concentration'}, hd, 1);
abundance.Statistic = abundance.Measure;
abundance.Replicate = ones(hd, 1);
abundance.Observation = cell(hd, 1); abundance.Observation(:) = {''};

abundance = abundance(:,allColumns);

Suaria_2020.abundance = abundance;
Suaria_2020.categories = mp_cat;
Suaria_2020.categories_macro = macro;
Suaria_2020.colour = mp_cat;
Suaria_2020.size = mp_size;

Data.Suaria_2020 = Suaria_2020;



%% Cunningham 2022

% Presence/absence data => no numeric values

source = 'Cunningham_2022';
filepath = fullfile(dataDirectory, source);

filename = 'raw data table S1.csv';
raw = readtable(fullfile(filepath, filename));

filename = 'air samples table S2.csv';
air = readtable(fullfile(filepath, filename));

filename = 'water samples table S3.csv';
water = readtable(fullfile(filepath, filename));

filename = 'sediment samples table S4.csv';
sediment = readtable(fullfile(filepath, filename));

filename = 'sea ice samples table S5.csv';
ice = readtable(fullfile(filepath, filename));

air.SampleType = cell(height(air), 1); air.SampleType(:) = {'air'};
water.SampleType = cell(height(water), 1); water.SampleType(:) = {'seawater'};
sediment.SampleType = cell(height(sediment), 1); sediment.SampleType(:) = {'sediment'};
ice.SampleType = cell(height(ice), 1); ice.SampleType(:) = {'ice'};

air.SampleID = air.FilterID;
water.SampleID = water.StationID;
sediment.SampleID = sediment.StationID;
ice.SampleID = ice.CoreID;

air.Date = [air.DeploymentDate, air.CollectionDate];
water.Date = [water.Date, NaT(height(water),1)];
sediment.Date = repmat(reshape(...
    datetime([2018 12 8; 2019 3 14], 'Format', 'yyyy-MM-dd'), 1, []), ...
    height(sediment), 1);
ice.Date = [ice.Date, NaT(height(ice), 1)];

air.Depth = cell(height(air), 1); air.Depth(:) = {''};
water.Depth = cell(height(water), 1); water.Depth(:) = {'subsurface'};
sediment.Depth = cellstr(num2str(sediment.Depth_m_));
ice.Depth = cell(height(ice), 1); ice.Depth(:) = {''};

air.Longitude_start = -air.Longitude_W_Start;
air.Longitude_end = -air.Longitude_W_End;
air.Latitude_start = -air.Latitude_S_Start;
air.Latitude_end = -air.Latitude_S_End;
air.Longitude = 0.5 .* (air.Longitude_start + air.Longitude_end);
air.Latitude = 0.5 .* (air.Latitude_start + air.Latitude_end);

water.Longitude = -water.Longitude_W_;
water.Latitude = -water.Latitude_S_;
water.Longitude_start = nan(height(water), 1); water.Longitude_end = nan(height(water), 1);
water.Latitude_start = nan(height(water), 1); water.Latitude_end = nan(height(water), 1);

sediment.Longitude = -sediment.Longitude_W_;
sediment.Latitude = -sediment.Latitude_S_;
sediment.Longitude_start = nan(height(sediment), 1); sediment.Longitude_end = nan(height(sediment), 1);
sediment.Latitude_start = nan(height(sediment), 1); sediment.Latitude_end = nan(height(sediment), 1);

ice.Longitude = -ice.Longitude_W_;
ice.Latitude = -ice.Latitude_S_;
ice.Longitude_start = nan(height(ice), 1); ice.Longitude_end = nan(height(ice), 1);
ice.Latitude_start = nan(height(ice), 1); ice.Latitude_end = nan(height(ice), 1);

water.Site = cell(height(water), 1); water.Site(:) = {''};
air.Site = cell(height(air), 1); air.Site(:) = {''};
ice.Site = cell(height(ice), 1); ice.Site(:) = {''};
sediment.Site = sediment.Location;

% combine data sets
getVars = {'SampleType', 'Site', 'SampleID', 'Date', 'Longitude', 'Latitude', ...
    'Longitude_start', 'Latitude_start', 'Longitude_end', 'Latitude_end', ...
    'Depth'};
dat = [water(:,getVars); air(:,getVars); ice(:,getVars); sediment(:,getVars)];

hd = height(dat);
dat.Source = repmat({[strrep(source, '_', ' ('), ')']}, hd, 1);
dat.SampleGear = cell(hd, 1); dat.SampleGear(:) = {''};
dat.LitterIDMethod = cell(hd, 1); dat.LitterIDMethod(:) = {''};
dat.LitterCategory = repmat({'plastic'}, hd, 1);
dat.LitterScale = cell(hd, 1); dat.LitterScale(:) = {''};
dat.PlasticForm = cell(hd, 1); dat.PlasticForm(:) = {'fibre'};
dat.PlasticSize = cell(hd, 1); dat.PlasticSize(:) = {''};
dat.SampleAtStation = isnan(dat.Longitude_start);
dat.SiteCategory = cell(hd, 1); dat.SiteCategory(:) = {''};
dat.Variable = cell(hd, 1); dat.Variable(:) = {'presence/absence'};
dat.Statistic = cell(hd, 1); dat.Statistic(:) = {'none'};
dat.Replicate = ones(hd, 1);
dat.Value = nan(hd, 1);
dat.Unit = cell(hd, 1); dat.Unit(:) = {'none'};
dat.Observation = cell(hd, 1); dat.Observation(:) = {'present'};

dat = dat(:,allColumns);

Cunningham_2022.abundance = dat;
Cunningham_2022.raw = raw;

Data.Cunningham_2022 = Cunningham_2022;




%% Combine abundance tables
sources = fieldnames(Data);
ns = length(sources);
nr = nan(1, ns);
nc = nan(1, ns);
for i = 1:ns
    source = sources{i};
    d = Data.(source).abundance;
    nr(i) = height(d);
    nc(i) = width(d);
end
nc = unique(nc);
if length(nc) > 1, warning([num2str(length(nc)-1) ' data sets have unexpected number of columns!']); end

DATA = table('Size', [sum(nr), nc], 'VariableTypes', struct2cell(columnTypes));
DATA.Properties.VariableNames = allColumns;
DATA.Date = [DATA.Date, DATA.Date];

cnr = [0 cumsum(nr)];
for i = 1:ns
    source = sources{i};
    d = Data.(source).abundance;
    j = cnr(i)+1:cnr(i+1);
    DATA(j,:) = d;
end







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

