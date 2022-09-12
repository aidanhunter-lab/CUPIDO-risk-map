%% Generate risk map

% Build up a risk map gradually layer by layer. Focus on a single year...


%% Preamble
% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('run_risk_map.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))

%% Model domain

% Specify the modelled area and horizontal & vertical resolutions
domain = map_domain();



%% Combine the next few code sections into a forcing data function...
%% Load GLORYS physical model -- monthly means

% Directory for model files
path_physicalModel = fullfile(baseDirectory, 'data', 'physical_models', ...
    'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern Ocean');
% Choose a single season
season = 2019;
% and months within that season
% months = [10:12, 1:3];
months = 1:3;
% Load the seawater density values.
% Model extraction & density calculations must be done already.
phys = loadWaterDensity(path_physicalModel, season, months);


% Match the data resolution to the map grid
% Interpolate over depths.
phys.density = permute(phys.density, [3 1 2 4]);
phys.density = interp1(phys.depth, phys.density, domain.depth);
phys.density = permute(phys.density, [2 3 1 4]);
phys.depth = domain.depth;
% Average over horizontal grid cells.
ncells = domain.nlon * domain.nlat;
ind = nan(length(phys.lon), length(phys.lat));
m = nan(ncells,1);
for j = 1:domain.nlat
    indj = domain.latgrid(j) <= phys.lat & phys.lat < domain.latgrid(j+1);
    for i = 1:domain.nlon
        indi = domain.longrid(i) <= phys.lon & phys.lon < domain.longrid(i+1);
        ind_ = indi & indj';
        n = (j - 1) * domain.nlon + i;
        m(n) = sum(ind_(:));
        ind(ind_) = n;
    end
end
nd = length(phys.depth);
nt = length(phys.time);
densityAv = nan(max(m), nd, nt, ncells);
for i = 1:ncells
    j = ind == i;
    j = repmat(j, [1 1 nd nt]);
    x = phys.density(j);
    x = reshape(x, m(i), nd, nt);
    densityAv(1:m(i), :, :, i) = x;
%     disp([num2str(round(100 * i / ncells, 2, 'significant')) '%'])
end

densityAv = mean(densityAv, 'omitnan');
densityAv = permute(densityAv, [4, 2, 3, 1]);
densityAv = reshape(densityAv, domain.nlon, domain.nlat, nd, nt);

phys.lon = domain.lon;
phys.lat = domain.lat;
phys.density = densityAv;

% Use the physical model data to determine which grid cells are land and
% which are below the seafloor.
domain.aboveSeafloor = ~isnan(phys.density(:,:,:,1)); % index all relevant lon-lat-depth cells
domain.isOcean = any(domain.aboveSeafloor,3); % index all relevant lon-lat cells


%% Load krill data

filename = 'krillbase_data.csv';
krill = readtable(filename);

% Map bounding vetices
mv = [domain.lon_range(1), domain.lon_range(1), domain.lon_range(2), ...
    domain.lon_range(2), domain.lon_range(1); ...
    domain.lat_range(1), domain.lat_range(2), domain.lat_range(2), ...
    domain.lat_range(1), domain.lat_range(1)];
% Find data within mapped region
inmap = inpolygon(krill.LONGITUDE, krill.LATITUDE, mv(1,:), mv(2,:));
% Omit data outside mapped region
krill = krill(inmap,:);

% Split date into years, months and days
[yr, mo, day] = datevec(krill.DATE);
krill.year = yr;
krill.month = mo;
krill.day = day;

% Choose a year range -- roughly corresponding to availability of other
% data, or simply use all available data
% ymin = 1997;
ymin = min(krill.year);
ymax = max(krill.year);
krill = krill(krill.year >= ymin & krill.year <= ymax,:);

% Use day-time and/or night-time samples
dayornight = 'both'; % select 'day', 'night', or 'both'
krill = krill(ismember(krill.DAY_NIGHT, {'day', 'night'}),:);
switch dayornight
    case {'day', 'night'}
        krill = krill(strcmp(krill.DAY_NIGHT, dayornight),:);
end

% Variable to plot
Var = 'STANDARDISED_KRILL_UNDER_1M2';

% Omit measuremtns zeros?
omitZeros = false;
switch omitZeros, case true
    krill = krill(krill.(Var) > 0,:);
end

krill_s = table2struct(krill, 'ToScalar', true);
fields = fieldnames(krill_s);
keepFields = {'LATITUDE', 'LONGITUDE', 'SEASON', 'TOP_SAMPLING_DEPTH', ...
    'BOTTOM_SAMPLING_DEPTH', 'STANDARDISED_KRILL_UNDER_1M2', 'year', ... 
    'month', 'day'};
krill_s = rmfield(krill_s, fields(~ismember(fieldnames(krill_s), keepFields)));

% Create struct for krill measurements -- similar to phys
zoo = rmfield(phys, 'density');
zoo.krill = zeros(domain.nlon, domain.nlat, length(domain.depth), length(phys.time));

% Smooth the measurements into a regular grid to match map domain.
% See Atkinson et al. 2008.
datgrid = nan(domain.nlon, domain.nlat, 1, length(phys.time));
for i = 1:domain.nlon
    indi = domain.longrid(i) < krill.LONGITUDE & krill.LONGITUDE <= domain.longrid(i+1);
    for j = 1:domain.nlat
        indj = indi & ...
            domain.latgrid(j) < krill.LATITUDE & krill.LATITUDE <= domain.latgrid(j+1);
        for k = 1:length(phys.time)
            tt = double(phys.time(k));
            [~,m] = datevec(tt);
            indk = indj & krill.month == m;
            if ~any(indk), continue; end
            dat = krill(indk,:);
            v = mean(dat.STANDARDISED_KRILL_UNDER_1M2, 'omitnan');
            datgrid(i,j,1,k) = v;
        end
    end
end

% Distribute krill between modelled depth layers.
% For now just put them all into the surface layer.
zoo.krill(:,:,1,:) = datgrid;

% Change units to number/m^3
zoo.krill =  zoo.krill .* ... 
    (0.5 * (domain.area(:,:,1:end-1) + domain.area(:,:,2:end))) ./ domain.volume;



%% Load phytoplankton data


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


















