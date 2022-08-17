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

% Number of lat-lon grid cells
nlon = length(waterDensity.lon);
nlat = length(waterDensity.lat);

%% Set up depth layers
% The GLORYS model output contains 50 depth layers going to nearly 6000 m,
% where the surface layers have fine resolution and the deep layers have
% resolution in the hundreds of meters.
% Let's choose our own depth scale with fewer layers. We are interested in
% the biologically active surface ocean; below about 500 m sinking faecal
% pellets may be considered exported to depth/free from biological
% recycling. For numerical stability, our model depth layers should not be
% too narrow... there's a trade-off that's related to the faecal pellet
% sinking speeds... the narrower the depth layers the smaller the model
% time steps...

% Particles sinking below exportDepth are assumed to be beyond
% remineralising microbial action
exportDepth = 500;

% Specify values for 'ndepths', 'w', and 'wfactor' that result in a
% sensible depth layer grid.
ndepths = 14; % choose the number of modelled depth layers.
w = 10; % width of uppermost depth layer
wfactor = 1.5; % depth layer widths increase with depth incrementally by 100*(wfactor-1)% 
w(2:ndepths,:) = wfactor .^ (1:ndepths-1) * w; % depth layer widths
wi = cumsum([0; w]);                           % depth layer edges
wm = 0.5 * (wi(1:end-1) + wi(2:end));          % depth layer midpoints

% The GLORYS model output corresponds to depth layer midpoints. Interpolate
% the model to match our selected depth layers.
waterDensity.density = permute(waterDensity.density, [3 1 2 4]); % shift depth to first dimension
% CHECK THE INTERP OPTIONS TO MAKE SURE THE BOUNDARIES ARE HANDLED WELL...
waterDensity.density = interp1(waterDensity.depth, waterDensity.density, wm);
waterDensity.density = permute(waterDensity.density, [2 3 1 4]); % restore dimension order
waterDensity.depth = wm;

% Index the ocean grid cells -- for efficiency, other cells may be omitted
% from calculations.
aboveSeafloor = ~isnan(waterDensity.density(:,:,:,1)); % index all relevant lon-lat-depth cells
isOceanCell = any(aboveSeafloor, 3); % index all relevant lon-lat cells

%% Calculate sinking speeds
% The sinking speed equations are written in cgs units.
shape = 'ellipsoid';  % faecal pellet shape
% rho_p = 1060 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3]
rho_p = 1300 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3]
mu = 0.00189 * 1e1;  % seawater viscosity, g / cm / s. [1 N s / m^2 = 1 Pa s = 10 g / cm / s]
g = 9.81 * 1e2;      % acceleration due to gravity, cm / s^2. [1 m / s^2 = 10^2 cm / s^2]

% Faecal pellet dimensions, cm.
% For now, let's assume that faecal pellets are either cylindrical OR
% ellipsoidal. We may introduce a mixture of shapes later...
% The values have been crudely approximated from Komar et al. 1981.
switch shape
    case 'cylinder'
        L = 0.09847; % length
        D = 0.01969; % diameter
        Ds = []; Di = []; Dl = [];
        V = 0.25 * pi * D ^ 2 * L; % volume
    case 'ellipsoid'
        Dl = 0.09847; % largest principal axis (diameter)
        Di = 0.03411; % intermediate
        Ds = 0.01706; % smallest
        L = []; D = [];
        V = pi/6 * Dl * Di * Ds; % volume
end

% Find faecal pellet sinking speed, cm / s
[v, Re] = sinking_speed(shape, 1e-3 * waterDensity.density, rho_p, mu, g, ...
    'L', L, 'D', D, 'Dl', Dl, 'Di', Di, 'Ds', Ds, ...
    'useMeans_Re', true, 'returnMax_Re', true);



%% Calculate faecal pellet production rates


















