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
dat = loadWaterDensity(path_physicalModel, season, months);

% Number of lat-lon grid cells
nlon = length(dat.lon);
nlat = length(dat.lat);

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
dat.density = permute(dat.density, [3 1 2 4]); % shift depth to first dimension
dat.density = interp1(dat.depth, dat.density, wm);
dat.density = permute(dat.density, [2 3 1 4]); % restore dimension order
dat.depth = wm;

% Index the ocean grid cells -- for efficiency, other cells may be omitted
% from calculations.
aboveSeafloor = ~isnan(dat.density(:,:,:,1)); % index all relevant lon-lat-depth cells
isOceanCell = any(aboveSeafloor, 3); % index all relevant lon-lat cells

%% Calculate sinking speeds
% The sinking speed equations are written in cgs units.
% Specify parameter values
pars.mu = 0.00189 * 1e1;  % seawater viscosity, g / cm / s. [1 N s / m^2 = 1 Pa s = 10 g / cm / s]
pars.g = 9.81 * 1e2;      % acceleration due to gravity, cm / s^2. [1 m / s^2 = 10^2 cm / s^2]
% See Atkinson et al., 2012, for krill faecal pellet properties
pars.rho_p = 1116 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3]
% pars.rho_p = 1220 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3] (value from Komar 1980)

% Faecal pellet dimensions (cm) crudely approximated from Komar et al 1981.
% Range of euphausiid faecal pellet volume measured by Komar is approx
% 1e-5 -> 5e-5. Antarctic krill likely near the upper side of this range.
% This is corroborated by Schmidt et al 2012, who use a standard pellet
% volume of 0.05 mm^3 = 5e-5 cm^3.
% pars.V = 5e-5;

% Atkinson et al 2012 provide more detail on krill faecal pellet
% dimensions, giving a larger median volume.
% pars.V = 1e-3 * 0.0754; % cm^3

pars.shape = 'cylinder';  % faecal pellet shape (may be cylinder or ellipsoid)
% It seems that the Atkinson data refer to cylindrical pellets
pars.L = 2667 * 1e-4; % median length (cm)
pars.D = 178 * 1e-4; % median width (cm)

% Faecal pellet dimensions, cm.
% For now, let's assume that faecal pellets are either cylindrical OR
% ellipsoidal. We may introduce a mixture of shapes later...
switch pars.shape
    case 'cylinder'
        pars.Ds = []; pars.Di = []; pars.Dl = [];
        pars.V = 0.25 * pi * pars.D ^ 2 * pars.L;
%         x = 15; % assume that length is x times diameter (L = xD)
%         pars.D = ((4 * pars.V) / (pi * x)) ^ (1 / 3); % diameter
%         pars.L = x * pars.D; % length
%         pars.Ds = []; pars.Di = []; pars.Dl = [];
    case 'ellipsoid'
        pars.Ds = pars.D; pars.Di = pars.D;
        pars.Dl = pars.L;
        pars.V = 1 / 6 * pi * pars.Ds * pars.Di * pars.Dl;
%         % assume that the two smallest diameters are equal (Di = Ds), and
%         % that the largest diameter is x times the smallest (Dl = xDi = xDs)
%         x = 5;
%         pars.Ds = ((6 * pars.V) / (pi * x)) ^ (1 / 3); % smallest principal axis (diameter)
%         pars.Di = pars.Ds; % intermediate
%         pars.Dl = x * pars.Ds; % largest
%         pars.L = []; pars.D = [];
end

% Find faecal pellet sinking speed, cm / s
[v, Re] = sinking_speed(pars.shape, 1e-3 * dat.density, pars.rho_p, pars.mu, pars.g, ...
    'L', pars.L, 'D', pars.D, 'Dl', pars.Dl, 'Di', pars.Di, 'Ds', pars.Ds, ...
    'useMeans_Re', true, 'returnMax_Re', true);



%% Include plastics

% Plastic dimensions. These may be comparable to the faecal pellet
% dimensions. Look up a paper that shows plastic particles inside faecal
% pellet and get approx dimensions from there.
% The plastics that are excreted may be smaller than those that are eaten
% because they are ground down to smaller sizes. An ingested piece of
% plastic may not be excreted all at once, it be excreted gradually as
% smaller pieces.

Vp = 229972525; % mu m^3
1e-12 * Vp




%% Calculate faecal pellet production rates


















