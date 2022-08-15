%% Generate risk map

% Build up a risk map gradually layer by layer. Focus on a single year...


%% Preamble
% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('risk_map.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))
path_physicalModel = fullfile(baseDirectory, 'data', 'physical_models', ... 
    'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern Ocean');



%% Load GLORYS physical model -- this should be wrapped into a function
season = 2019;
months = [10:12, 1:3];
waterDensity = loadWaterDensity(path_physicalModel, season, months);



%% Calculate sinking speeds
shape = 'cylinder';
rho_p = 1100; % faecal pellet density
mu = 0.0189; % g / cm / s -- seawater viscosity [1 N s / m^2 = 1 Pa s = 10 g / cm / s]
g = 981; % cm / s^2;


% I need to figure out the units here...
v = sinking_speed(shape, waterDensity.density, rho_p, mu, g, ...
    'L', 0.01, 'D', 0.005);



%% Calculate faecal pellet production rates


















