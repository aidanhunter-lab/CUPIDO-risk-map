%% Load and organise forcing data
function [forc, domain] = prepare_forcing(baseDirectory, domain, varargin)


%% Set default values for optional arguments

extractVarargin(varargin)

% Directory for model files
if ~exist('path_physicalModel', 'var')
    path_physicalModel = fullfile(baseDirectory, 'data', 'physical_models', ...
        'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern Ocean');
end
% Choose a single season
if ~exist('season', 'var')
    % data for selected season must be stored in directory path_physicalModel
    season = 2019;
end
% and months within that season
months = domain.month';
% if ~exist('months', 'var')
%     months = 1:3; % choose Jan-Mar to match available krill data
% %     months = [10:12, 1:3];
% end

% Years of krillbase data to retain -- 2-element vector (min,max)
if ~exist('year_range_krillbase', 'var')
    year_range_krillbase = [-inf,inf]; % All years by default
end

% Use day-time and/or night-time samples from krillabse
if ~exist('dayornight_krillbase', 'var')
    dayornight_krillbase = 'both'; % select 'day', 'night', or 'both'
end
if ~ismember(dayornight_krillbase, {'both','day','night'})
    error("Input argument dayornight_krillbase must be 'both', 'day', or 'night'.")
end

% Krillbase measurement variable to use as abundance value
if ~exist('var_krillbase', 'var')
    var_krillbase = 'STANDARDISED_KRILL_UNDER_1M2';
end


% SeaWiFS chlorophyll measurement variables to use as abundance value(s)
if ~exist('vars_chl', 'var')
    vars_chl = 'Tchl'; % total chlorophyll as default -- returning multiple measurement variables may be possible
end


%% GLORYS physical model -- monthly means

% Load the seawater density values.
% Model extraction & density calculations must be done already.
phys = loadWaterDensity(path_physicalModel, season, months);

% Match the data resolution to the map grid
[~,phys.month] = datevec(double(phys.time));

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

nd = domain.ndepth;
nt = domain.nmonth;
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
% Include index matrices in the domain struct
isLand = all(isnan(phys.density(:,:,:,1)), 3);
isWaterColumn = ~isnan(phys.density(:,:,:,1));
isSurfaceLayer = cat(3, true(domain.nlon, domain.nlat, 1), ... 
    false(domain.nlon, domain.nlat, domain.ndepth-1)) .* ~isLand;
isBottomLayer = cat(3, false(domain.nlon, domain.nlat), ...
    logical(-diff(isWaterColumn, 1, 3)));

domain.isLand = isLand;
domain.isWaterColumn = isWaterColumn;
domain.isSurfaceLayer = isSurfaceLayer;
domain.isBottomLayer = isBottomLayer;


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

% Filter by year -- by default all data is retained
krill = krill(krill.year >= year_range_krillbase(1) & ... 
    krill.year <= year_range_krillbase(2),:);

% Use day-time and/or night-time samples
krill = krill(ismember(krill.DAY_NIGHT, {'day', 'night'}),:); % omit samples where day or night is unknown
switch dayornight_krillbase
    case {'day', 'night'}
        krill = krill(strcmp(krill.DAY_NIGHT, dayornight_krillbase),:);
end

% Omit measuremtns zeros?
omitZeros = false;
switch omitZeros, case true
    krill = krill(krill.(var_krillbase) > 0,:);
end

% Create struct for krill measurements -- similar to phys
zoo = rmfield(phys, 'density');
zoo.krill = zeros(domain.nlon, domain.nlat, domain.ndepth, domain.nmonth);

% Smooth the measurements into a regular grid to match map domain.
% See Atkinson et al. 2008.
datgrid = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
for i = 1:domain.nlon
    indi = domain.longrid(i) < krill.LONGITUDE & krill.LONGITUDE <= domain.longrid(i+1);
    for j = 1:domain.nlat
        indj = indi & ...
            domain.latgrid(j) < krill.LATITUDE & krill.LATITUDE <= domain.latgrid(j+1);
        for k = 1:domain.nmonth
            m = domain.month(k);
            indk = indj & krill.month == m;
            if ~any(indk), continue; end
            dat = krill(indk,:);
            v = mean(dat.(var_krillbase), 'omitnan');
            datgrid(i,j,1,k) = v;
        end
    end
end

% Index grid cells lacking data
noKrillData = isnan(datgrid);
% and set all these as NaN
zoo.krill(repmat(noKrillData, [1, 1, domain.ndepth, 1])) = nan;
% Distribute observed krill between modelled depth layers.
% For now just put them all into the surface layer.
zoo.krill(:,:,1,:) = datgrid;

% Change units to number/m^3
zoo.krill =  zoo.krill .* ... 
    (0.5 * (domain.area(:,:,1:end-1) + domain.area(:,:,2:end))) ./ domain.volume;


%% Phytoplankton data

filename = 'SeaWiFS_Phytoplankton-Size-Class-1997-2007-Southern-Ocean.mat';
load(filename, 'Dat') % data was stored as variable 'Dat'
% Convert -999 to nan
fields = fieldnames(Dat);
for i = 1:length(fields)
    x = Dat.(fields{i});
    x(x == -999) = nan;
    Dat.(fields{i}) = x;
end

% Omit unused variables
fields = fieldnames(Dat);
keepVars = [{'month','latitude','longitude'}, vars_chl];
Dat = rmfield(Dat, fields(~ismember(fields, keepVars)));

% Filter the data by month
ind = ismember(Dat.month, months);
Dat = structfun(@(z) z(ind), Dat, ...
    'UniformOutput', false);

% Find data within mapped region
inmap = inpolygon(Dat.longitude, Dat.latitude, mv(1,:), mv(2,:));
% Omit data outside mapped region
Dat = structfun(@(z) z(inmap), Dat, 'UniformOutput', false);

% Create struct for phytoplankton prey measurements -- similar to phys
prey = rmfield(phys, 'density');
prey.chl = zeros(domain.nlon, domain.nlat, domain.ndepth, domain.nmonth);

% Smooth the measurements into a regular grid to match map domain.
datgrid = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
for i = 1:domain.nlon
    indi = domain.longrid(i) < Dat.longitude & Dat.longitude <= domain.longrid(i+1);
    for j = 1:domain.nlat
        indj = indi & ...
            domain.latgrid(j) < Dat.latitude & Dat.latitude <= domain.latgrid(j+1);
        for k = 1:domain.nmonth
            m = domain.month(k);
            indk = indj & Dat.month == m;
            if ~any(indk), continue; end
            dat = structfun(@(z) z(indk), Dat, 'UniformOutput', false);
            v = mean(dat.(vars_chl), 'omitnan'); % this is not robust to multiple useVars values
            datgrid(i,j,1,k) = v;
%             disp(((i-1)*domain.nlat*length(prey.time) + (j-1)*length(prey.time) + k) / domain.nlon / domain.nlat / length(prey.time))
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



%% Plastics


%% Output struct
forc.lon = domain.lon;
forc.lat = domain.lat;
forc.depth = domain.depth;
forc.month = domain.month;
forc.time = double(phys.time);

forc.density_seawater = phys.density;
forc.chl_total = prey.chl;
forc.krill = zoo.krill;
forc.noKrillData = noKrillData;

forcSize = [domain.nlon, domain.nlat, domain.ndepth, domain.nmonth];

correctSize = [all(size(forc.density_seawater) == forcSize), ...
    all(size(forc.chl_total) == forcSize), ... 
    all(size(forc.krill) == forcSize)];

if ~all(correctSize)
    warning('Some forcing data have incorrect dimension! Something is not right...')
end

