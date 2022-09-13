%% Load and organise forcing data
function [forc, domain] = prepare_forcing(baseDirectory, domain, varargin)

extractVarargin(varargin)

% Directory for model files
if ~exist('path_physicalModel', 'var')
    path_physicalModel = fullfile(baseDirectory, 'data', 'physical_models', ...
        'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern Ocean');
end
% Choose a single season
if ~exist('season', 'var')
    season = 2019;
end
% and months within that season
if ~exist('months', 'var')
    months = 1:3; % choose Jan-Mar to match available krill data
%     months = [10:12, 1:3];
end


%% GLORYS physical model -- monthly means

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

forc.phys = phys;


%% Krill data

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

forc.zoo = zoo;


%% Phytoplankton data

filename = 'SeaWiFS_Phytoplankton-Size-Class-1997-2007-Southern-Ocean.mat';
load(filename, 'Dat') % data was stored as variable 'Dat'
% disp(Dat)
% Convert -999 to nan
fields = fieldnames(Dat);
for i = 1:length(fields)
    x = Dat.(fields{i});
    x(x == -999) = nan;
    Dat.(fields{i}) = x; clear x
end

% Omit unused variables
useVars = 'Tchl';
fields = fieldnames(Dat);
keepVars = [{'month','latitude','longitude'}, useVars];
Dat = rmfield(Dat, fields(~ismember(fields, keepVars)));

% Filter the data by month
ind = ismember(Dat.month, months);
Dat = structfun(@(z) z(ind), Dat, ...
    'UniformOutput', false);

% Exclude data outside map boundaries
mv = [domain.lon_range(1), domain.lon_range(1), domain.lon_range(2), ...
    domain.lon_range(2), domain.lon_range(1); ...
    domain.lat_range(1), domain.lat_range(2), domain.lat_range(2), ...
    domain.lat_range(1), domain.lat_range(1)];
% Find data within mapped region
inmap = inpolygon(Dat.longitude, Dat.latitude, mv(1,:), mv(2,:));
% Omit data outside mapped region
Dat = structfun(@(z) z(inmap), Dat, 'UniformOutput', false);

% Create struct for phytoplankton prey measurements -- similar to phys
prey = rmfield(phys, 'density');
prey.chl = zeros(domain.nlon, domain.nlat, length(domain.depth), length(phys.time));

% Smooth the measurements into a regular grid to match map domain.
datgrid = nan(domain.nlon, domain.nlat, 1, length(prey.time));
for i = 1:domain.nlon
    indi = domain.longrid(i) < Dat.longitude & Dat.longitude <= domain.longrid(i+1);
    for j = 1:domain.nlat
        indj = indi & ...
            domain.latgrid(j) < Dat.latitude & Dat.latitude <= domain.latgrid(j+1);
        for k = 1:length(prey.time)
            tt = double(prey.time(k));
            [~,m] = datevec(tt);
            indk = indj & Dat.month == m;
            if ~any(indk), continue; end
            dat = structfun(@(z) z(indk), Dat, 'UniformOutput', false);
            v = mean(dat.(useVars), 'omitnan'); % this is not robust to multiple useVars values
            datgrid(i,j,1,k) = v;
            disp(((i-1)*domain.nlat*length(prey.time) + (j-1)*length(prey.time) + k) / domain.nlon / domain.nlat / length(prey.time))
        end
    end
end


% Distribute prey between modelled depth layers.
% For now just put them all into the surface layer.
prey.chl(:,:,1,:) = datgrid;

% FIGURE OUT THE UNITS!!

% % Change units to number/m^3
% zoo.krill =  zoo.krill .* ... 
%     (0.5 * (domain.area(:,:,1:end-1) + domain.area(:,:,2:end))) ./ domain.volume;

forc.prey = prey;

