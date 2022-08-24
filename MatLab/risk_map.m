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
% Load parameter values
pars = initialise_parameters();

% Find faecal pellet sinking speed, cm / s
[v, Re] = sinking_speed(pars.shape, 1e-3 * dat.density, pars.rho_p, pars.mu, pars.g, ...
    'L', pars.L, 'D', pars.D, 'Dl', pars.Dl, 'Di', pars.Di, 'Ds', pars.Ds, ...
    'useMeans_Re', false, 'returnMax_Re', true);

% PERHAPS THE SINKING SPEED EQUATIONS COULD BE MODIFIED TO ACCOUNT FOR
% VOLUME REDUCTIONS RESULTING FROM REMINERALISATION. I DON'T WANT TO MODEL
% FAECAL PELLETS OF VARYING SIZE (THEY WILL BE CONSTANT VOLUME), HOWEVER, I
% DO WANT TO MODEL REMINERALISATION WHICH IS A PROCESS THAT REDUCES SIZE
% AND THEREFORE REDUCES SINK SPEED. REMINERALISATION WILL BE MODELLED AS
% REDUCTION IN TOTAL CARBON (REDUCING THE NUMBER OF FAECAL PELLETS, BUT NOT
% THEIR SIZE). I COULD TRY TO MODIFY THE SINK SPEED EQUATION TO ACCOUNT FOR
% THIS BY REDUCING SINK SPEED WITHIN DEEPER LAYERS WHERE PELLETS HAVE SPENT
% MORE TIME BEING DEGRADED. THIS WILL TAKE A BIT OF THOUGHT, BUT IS WORTH
% KEEPING IN MIND... I'M NOT SURE HOW LARGE/IMPORTANT THE EFFECT WILL BE...


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


















