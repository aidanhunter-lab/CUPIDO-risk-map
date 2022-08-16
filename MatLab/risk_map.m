%% Generate risk map

% Build up a risk map gradually layer by layer. Focus on a single year...


%% Preamble
% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('risk_map.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))

%% Load GLORYS physical model -- monthly means

% Directory for model files
path_physicalModel = fullfile(baseDirectory, 'data', 'physical_models', ...
    'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern Ocean');
% Choose a single season
season = 2019;
% and months within that season
months = [10:12, 1:3];
% Load the seawater density values.
% Model extraction & density calculations must be done already.
waterDensity = loadWaterDensity(path_physicalModel, season, months);


%% Calculate sinking speeds
% The sinking speed equations are written in cgs units.
shape = 'cylinder';  % faecal pellet shape
rho_p = 1060 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3]
mu = 0.00189 * 1e1;  % seawater viscosity, g / cm / s. [1 N s / m^2 = 1 Pa s = 10 g / cm / s]
g = 9.81 * 1e2;      % acceleration due to gravity, cm / s^2. [1 m / s^2 = 10^2 cm / s^2]

% Faecal pellet dimensions, cm.
% For now, let's assume that faecal pellets are either cylindrical OR
% ellipsoidal. We may introduce a mixture of shapes later...
switch shape
    case 'cylinder'
        L = 0.25; % length
        D = 0.05; % cross-sectional diameter
        Ds = []; Di = []; Dl = [];
    case 'ellipsoid'
        Dl = 0.1; % largest principal axis (diameter)
        Di = 0.05; % intermediate
        Ds = 0.025; % smallest
        L = []; D = [];
end

% Find faecal pellet sinking speed, cm / s
v = sinking_speed(shape, 1e-3 * waterDensity.density, rho_p, mu, g, ...
    'L', L, 'D', D, 'Dl', Dl, 'Di', Di, 'Ds', Ds);



%% Calculate faecal pellet production rates


















