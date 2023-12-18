%% Prepare the plastic data

%% Directories
thisFile = which('prepare_plastic_data.m');
dirBase = fileparts(fileparts(thisFile));
dirData = fullfile(dirBase, 'data');
dirDataSource = fullfile(dirData, 'plastic samples', 'sources');
dirDataSave = fullfile(dirData, 'plastic samples', 'collated');
dirMatLab = fullfile(dirBase, 'MatLab');
addpath(genpath(dirData))
addpath(genpath(dirMatLab))

%% Load data and coerce into single table
[plastics, plastics_table] = load_plastic_conc_data(...
    'dataDirectory', dirDataSource);

%% Regularise & filter data & save
[abundance, sources, nsources] = filter_plastic_conc_data(...
    plastics, plastics_table, ...
    'saveData', true, 'outputFilename', 'plastic_quantity.csv', ...
    'FilterByLitterCategory', false, 'FilterByLitterScale', false, ...
    'regularisePlasticForm', false);
% if saveData=true then this collated data table is saved as .csv file
head(abundance)

