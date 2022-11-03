function pars = initialise_parameters(domain, varargin)
% Set default values for all model parameters.

% We may specify any parameter value in the function call as an optional
% argument name-value pair.
extractVarargin(varargin)

if ~exist('polymer', 'var')
    polymer = 'Low density Polyethylene'; % see parsFromFile (below) for options
end

%% Load parameters
if ~exist('loadParsFromFile', 'var')
    loadParsFromFile = true;
end
if ~exist('parsFileName', 'var')
    parsFileName = 'parameters.csv';
end
switch loadParsFromFile, case true
    if exist(parsFileName, 'file') == 2
        parsFromFile = readtable(parsFileName);
    else
        parsFromFile = [];
    end
    otherwise
        parsFromFile = [];
end
if isempty(parsFromFile)
    warning(['Parameters have not been loaded from file! Is ' parsFileName ' saved the data directory?'])
    loadParsFromFile = false;
end

%% Modelled time steps
if ~exist('dt_out', 'var')
    % time step (s) of returned results
    dt_out = 3600; % 60 minute outputs
end
if ~exist('dt_max', 'var')
    % maximum integration time step (s)
    dt_max = 1800; % 30 minute steps
end
% There is a stability constraint on integration time steps...
min_width = min(diff(domain.depthgrid)); % width (m) of narrowest depth layer
max_speed = 2000 / 24 / 60 / 60; % upper limit on faecal pellet sinking speed (m/s)
dt_max_ = min_width / max_speed;
if dt_max >= dt_max_
    % If time step is too large then reset to (hopefully!) stable value
    dt_max = dt_max_;
    dt_max = floor(dt_max / 60 / 5) * 5 * 60; % round down to nearest 5 minutes
end
pars.dt_out = dt_out;
pars.dt_max = dt_max;

%% Physical parameters
pars.mu = 0.00189 * 1e1;  % seawater viscosity, g / cm / s. [1 N s / m^2 = 1 Pa s = 10 g / cm / s]
pars.g = 9.81 * 1e2;      % acceleration due to gravity, cm / s^2. [1 m / s^2 = 10^2 cm / s^2]

% See Atkinson et al., 2012, for krill faecal pellet properties
pars.rho_f = 1116 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3]
% pars.rho_p = 1220 * 1e-3; % value from Komar 1980 

% For now, let's assume that faecal pellets are either cylindrical OR
% ellipsoidal. We may introduce a mixture of shapes later...
pars.shape = 'cylinder';  % faecal pellet shape (may be cylinder or ellipsoid)

% It seems that Atkinson's data refer to cylindrical pellets
pars.L = 2667 * 1e-4; % median length (cm)
pars.D = 178 * 1e-4; % median width (diameter) (cm)

switch pars.shape
    case 'cylinder'
        pars.Ds = []; pars.Di = []; pars.Dl = [];
        pars.V = 0.25 * pi * pars.D ^ 2 * pars.L;
    case 'ellipsoid'
        % assume that the two smallest diameters are equal
        pars.Ds = pars.D; pars.Di = pars.D;
        pars.Dl = pars.L;
        pars.V = 1 / 6 * pi * pars.Ds * pars.Di * pars.Dl;
end
% The faecal pellet volumes from the Atkinson data are at the uppper end of
% the range given for euphausiid pellets (approx 1e-5 -> 5e-5) in 
% Komar et al (1981). This makes sense as Antarctic krill are largest
% euphausiid. The volume is also close to the standard volume (5e-5) used
% by Schmidt et al (2012).

% Faecal pellet carbon content - measured by Clara Manno from sediment trap
% samples (original units: mg C / mm^3)
switch pars.shape
    case 'cylinder'
        pars.fpCm_summer = 0.03;   % mean g C / cm^3
        pars.fpCsd_summer = 0.006; % standard deviation
        pars.fpCm_winter = 0.018;
        pars.fpCsd_winter = 0.006;
    case 'ellipsoid'
        pars.fpCm_summer = 0.052;
        pars.fpCsd_summer = 0.005;
        pars.fpCm_winter = 0.034;
        pars.fpCsd_winter = 0.006;
    case 'sphere'
        pars.fpCm_summer = 0.035;
        pars.fpCsd_summer = 0.004;
        pars.fpCm_winter = 0.027;
        pars.fpCsd_winter = 0.008;
end


%% Krill size conversions
switch loadParsFromFile, case true
    pars.W_dry_a = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'W_dry_a'),:); % mg dry / mm
    pars.W_dry_b = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'W_dry_b'),:);
    pars.W_wet_a = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'W_wet_a'),:); % mg wet / mm
    pars.W_wet_b = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'W_wet_b'),:);
    pars.W_dry2c = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'W_dry2c') & ...
        strcmp(parsFromFile.Group, 'Summer'),:); % mg C / mg dry
end

%% Faecal pellet production rate
switch loadParsFromFile, case true
    pars.FPprod = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'FPprod')); % mm^3 / h / g dry mass
end


%% Plastic properties
switch loadParsFromFile, case true
    pars.Vp = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'V_m'),:); % mum^3 / particle
    pars.rho_p = parsFromFile.Value(strcmp(parsFromFile.Parameter, 'rho_m') & ...
        strcmp(parsFromFile.Group, polymer));
end

