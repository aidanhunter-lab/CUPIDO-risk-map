%% Load and organise forcing data
function [forc, domain, map_krill, map_krill_table, map_chl_table, map_chl_table_highRes] = ... 
    prepare_forcing(baseDirectory, domain, varargin)

extractVarargin(varargin)

%% Load forcing and domain from file?
if ~exist('loadFromFile', 'var')
    loadFromFile = false;
end
if ~exist('forcFile', 'var')
    forcFile = 'forcing.mat';
end
forcPath = fullfile(baseDirectory, 'MatLab', 'temp', forcFile);
if ~exist('domainFile', 'var')
    domainFile = 'domain.mat';
end
domainPath = fullfile(baseDirectory, 'MatLab', 'temp', domainFile);

switch loadFromFile, case true
    forc = load(forcPath).forc;
    domain = load(domainPath).domain;
    return
end

%% Set default values for optional arguments

% Save output struct in 'temp' directory to avoid recalling function?
if ~exist('saveOutput', 'var')
    saveOutput = true;
end

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

% Depth layers with krill
if ~exist('krillDepths', 'var')
    krillDepths = 1; % default is top depth layer only
end

% SeaWiFS chlorophyll measurement variables to use as abundance value(s)
if ~exist('vars_chl', 'var')
    vars_chl = 'Tchl'; % total chlorophyll as default -- returning multiple measurement variables may be possible
end

% Time forcing data set-up
if ~exist('timeSetUp', 'var')
    timeSetUp = true;
end
switch timeSetUp, case true
    tic; fprintf('\n\n'); disp(append("started: ", string(datetime('now'))))
end

if ~exist('makeHiResChlMap', 'var')
    makeHiResChlMap = false;
end

if ~exist('chlExtraResolution', 'var')
    chlExtraResolution = 3; % chl resolution is chlExtraResolution * krill resolution
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

% Density
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
A = nan(domain.nlon, domain.nlat, domain.ndepth, domain.nmonth);
zoo.krill = A;

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
% zoo.krill(repmat(noKrillData, [1, 1, domain.ndepth, 1])) = nan;

% Put all krill into top depth layer, later they may be distributed across depth
zoo.krill(:,:,1,:) = datgrid; % individuals/m^2
map_krill = squeeze(datgrid); % output for plotting map of krill abundance

% Change units to individuals/grid cell
zoo.krill =  zoo.krill .* ... 
    (0.5 * (domain.area(:,:,1:end-1) + domain.area(:,:,2:end))); % individuals/grid cell

% % Change units to number/m^3
% zoo.krill =  zoo.krill .* ... 
%     (0.5 * (domain.area(:,:,1:end-1) + domain.area(:,:,2:end))); % individuals/grid cell
% zoo.krill =  zoo.krill ./ domain.volume; % individuals/m^3


% Length
filename = 'LENGTH_DATA_20200201.csv';
krillLength = readtable(filename);

filename = 'HAUL_METADATA_20200201.csv';
meta = readtable(filename);

% filename = 'NET_CODES20200201.csv';
% netcodes = readtable(filename);
% 
% filename = 'NET_TRAJECTORY_20200201.csv';
% nettraj = readtable(filename);

krillLength = innerjoin(krillLength, meta, 'Keys', 'REF');

% Find data within mapped region
inmap = inpolygon(krillLength.LON, krillLength.LAT, mv(1,:), mv(2,:));
% Omit data outside mapped region
krillLength = krillLength(inmap,:);

% Filter by year -- by default all data is retained
krillLength = krillLength(krillLength.YEAR >= year_range_krillbase(1) & ... 
    krillLength.YEAR <= year_range_krillbase(2),:);

% Use day-time and/or night-time samples
krillLength = krillLength(ismember(krillLength.DAY_NIGHT, {'D', 'N'}),:); % omit samples where day or night is unknown
switch dayornight_krillbase
    case {'day', 'night'}
        krillLength = krillLength(strcmp(krillLength.DAY_NIGHT, dayornight_krillbase),:);
end

% Unit conversions

if exist('pars', 'var')
    pars = evalin('caller', 'pars');
    try
        % length to weight
        krillLength.WEIGHT_dry = pars.W_dry_a .* krillLength.LENGTH_mm_ .^ pars.W_dry_b;
        krillLength.WEIGHT_wet = pars.W_wet_a .* krillLength.LENGTH_mm_ .^ pars.W_wet_b;
        % carbon content from dry weight
        krillLength.CARBON = pars.W_dry2c .* krillLength.WEIGHT_dry;
    catch
        warning("Optional argument 'pars' is present but the allometric parameters are missing. See 'intialise_parameters.m'.")
    end
else
    filename = 'parameters.csv';
    try
        pars = readtable(filename);
    catch
        warning("Optional argument 'pars' was not specified, so parameter table (parameters.csv) has been loaded from file.")
    end
    % length to weight
    a = 'W_dry_a';
    b = 'W_dry_b';
    a = pars.Value(strcmp(pars.Parameter, a));
    b = pars.Value(strcmp(pars.Parameter, b));
    krillLength.WEIGHT_dry = a .* krillLength.LENGTH_mm_ .^ b;
    a = 'W_wet_a';
    b = 'W_wet_b';
    a = pars.Value(strcmp(pars.Parameter, a));
    b = pars.Value(strcmp(pars.Parameter, b));
    krillLength.WEIGHT_wet = a .* krillLength.LENGTH_mm_ .^ b;
    % carbon content from dry weight
    p = 'W_dry2c';
    s = 'Summer';
    p = pars.Value(strcmp(pars.Parameter, p) & strcmp(pars.Group, s));
    krillLength.CARBON = p .* krillLength.WEIGHT_dry;

    pars_ = pars(:,{'Parameter','Value'});
    clearvars pars
    for i = 1:height(pars_)
        pars.(pars_.Parameter{i}) = pars_.Value(i);
    end
end


% Estimate krill biomass from length data and allometric conversions
zoo.krillDryWeight = A;
zoo.krillWetWeight = A;
zoo.krillC = A;
minL = min(krillLength.LENGTH_mm_);
maxL = max(krillLength.LENGTH_mm_);
lens = minL:maxL;
dryWeights = pars.W_dry_a .* lens .^ pars.W_dry_b;
wetWeights = pars.W_wet_a .* lens .^ pars.W_wet_b;
CWeights = pars.W_dry2c .* dryWeights;
% wscale = 1e-3; % we want to return weight units of g not mg

for i = 1:domain.nlon
    indi = domain.longrid(i) < krillLength.LON & krillLength.LON <= domain.longrid(i+1);
    for j = 1:domain.nlat
        indj = indi & ...
            domain.latgrid(j) < krillLength.LAT & krillLength.LAT <= domain.latgrid(j+1);
        for k = 1:domain.nmonth
            m = domain.month(k);
            indk = indj & krillLength.MONTH == m;
            if ~any(indk), continue; end
            dat = krillLength(indk,:);
            l = dat.LENGTH_mm_;
            lsmooth = fitdist(l,'Kernel');
            y = pdf(lsmooth,lens);
            nz = zoo.krill(i,j,1,k); % individuals
            y = nz * y; % individuals at length
            dw = sum(y .* dryWeights); % dry weight g
            ww = sum(y .* wetWeights); % wet weight g
            cw = sum(y .* CWeights); % carbon weight g
            zoo.krillDryWeight(i,j,1,k) = dw;
            zoo.krillWetWeight(i,j,1,k) = ww;
            zoo.krillC(i,j,1,k) = cw;
        end
    end
end

% Calculate summary statistics for length and weight
zoo.krillLenMeanLog = A;
zoo.krillLenSDLog = A;
% dry weight
zoo.krillDryWtMeanLog = A;
zoo.krillDryWtSDLog = A;
% wet weight
zoo.krillWetWtMeanLog = A;
zoo.krillWetWtSDLog = A;
% carbon weight
zoo.krillCMeanLog = A;
zoo.krillCSDLog = A;

for i = 1:domain.nlon
    indi = domain.longrid(i) < krillLength.LON & krillLength.LON <= domain.longrid(i+1);
    for j = 1:domain.nlat
        indj = indi & ...
            domain.latgrid(j) < krillLength.LAT & krillLength.LAT <= domain.latgrid(j+1);
        for k = 1:domain.nmonth
            m = domain.month(k);
            indk = indj & krillLength.MONTH == m;
            if ~any(indk), continue; end
            dat = krillLength(indk,:);

            % Mean lengths weighted by event sample size
            ev = unique(dat.EVENT);
            evn = arrayfun(@(z) sum(dat.EVENT == z), ev); % number of measures per sample event
            % Summary stats on log scale
            meanLogLen_ = arrayfun(@(z) mean(log(dat.LENGTH_mm_(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogLen_ = sum(evn ./ sum(evn) .* meanLogLen_);
            sdLogLen_ = arrayfun(@(z) std(log(dat.LENGTH_mm_(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogLen_ = sum(evn ./ sum(evn) .* sdLogLen_);

            meanLogWdry_ = arrayfun(@(z) mean(log(dat.WEIGHT_dry(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogWdry_ = sum(evn ./ sum(evn) .* meanLogWdry_);
            sdLogWdry_ = arrayfun(@(z) std(log(dat.WEIGHT_dry(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogWdry_ = sum(evn ./ sum(evn) .* sdLogWdry_);

            meanLogWwet_ = arrayfun(@(z) mean(log(dat.WEIGHT_wet(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogWwet_ = sum(evn ./ sum(evn) .* meanLogWwet_);
            sdLogWwet_ = arrayfun(@(z) std(log(dat.WEIGHT_wet(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogWwet_ = sum(evn ./ sum(evn) .* sdLogWwet_);

            meanLogC_ = arrayfun(@(z) mean(log(dat.CARBON(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogC_ = sum(evn ./ sum(evn) .* meanLogC_);
            sdLogC_ = arrayfun(@(z) std(log(dat.CARBON(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogC_ = sum(evn ./ sum(evn) .* sdLogC_);

            zoo.krillLenMeanLog(i,j,:,k) = meanLogLen_;
            zoo.krillLenSDLog(i,j,:,k) = sdLogLen_;
            zoo.krillDryWtMeanLog(i,j,:,k) = meanLogWdry_;
            zoo.krillDryWtSDLog(i,j,:,k) = sdLogWdry_;
            zoo.krillWetWtMeanLog(i,j,:,k) = meanLogWwet_;
            zoo.krillWetWtSDLog(i,j,:,k) = sdLogWwet_;
            zoo.krillCMeanLog(i,j,:,k) = meanLogC_;
            zoo.krillCSDLog(i,j,:,k) = sdLogC_;
        end
    end
end


% Some density samples are not matched to length stats -- infill with
% nearest neighbours
gotL = ~isnan(zoo.krillLenMeanLog); % index data values
% gotL = ~isnan(zoo.krillDryWeight); % index data values
misL = ~isnan(zoo.krill) & ~gotL; % index data to infill
for j = 1:length(zoo.month)
    g = gotL(:,:,1,j);
    i = misL(:,:,1,j);
    li = find(i);
    [ri, ci] = find(i);
    fi = [ri, ci];
    [ri, ci] = find(g);
    gi = [ri, ci];
    ni = sum(i(:));
    for k = 1:ni
        nz = zoo.krill(:,:,1,j);
        % length
        ml = zoo.krillLenMeanLog(:,:,1,j);
        ms = zoo.krillLenSDLog(:,:,1,j);
        ii = fi(k,:);
        d = (sum(abs(gi - ii) .^ 2, 2)) .^ 0.5;
        md = find(d == min(d), 1);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillLenMeanLog(:,:,1,j) = ml;
        zoo.krillLenSDLog(:,:,1,j) = ms;
        % dry weight
        ml = zoo.krillDryWtMeanLog(:,:,1,j);
        ms = zoo.krillDryWtSDLog(:,:,1,j);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillDryWtMeanLog(:,:,1,j) = ml;
        zoo.krillDryWtSDLog(:,:,1,j) = ms;
        m = zoo.krillDryWeight(:,:,1,j);
        m(li(k)) = exp(ml(gi(md,1), gi(md,2))) * nz(gi(md,1), gi(md,2));
        zoo.krillDryWeight(:,:,1,j) = m;
        % wet weight
        ml = zoo.krillWetWtMeanLog(:,:,1,j);
        ms = zoo.krillWetWtSDLog(:,:,1,j);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillWetWtMeanLog(:,:,1,j) = ml;
        zoo.krillWetWtSDLog(:,:,1,j) = ms;
        m = zoo.krillWetWeight(:,:,1,j);
        m(li(k)) = exp(ml(gi(md,1), gi(md,2))) * nz(gi(md,1), gi(md,2));
        zoo.krillWetWeight(:,:,1,j) = m;
        % carbon weight
        ml = zoo.krillCMeanLog(:,:,1,j);
        ms = zoo.krillCSDLog(:,:,1,j);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillCMeanLog(:,:,1,j) = ml;
        zoo.krillCSDLog(:,:,1,j) = ms;
    end
end

% Distribute observed krill between modelled depth layers.
% Krill depths may be the top layer only, or multiple contiguous layers
% including the top layer. Krill are assumed to be evenly distributed
% between the selected depth layers -- some distribution function may be
% specified at later time...
% krillDepths = 1;
wn = length(krillDepths);
wd = false(domain.ndepth, 1); wd(krillDepths) = true;

zoo.krill(:,:,wd,:) = repmat(zoo.krill(:,:,1,:), [1 1 wn 1]) ./ wn;
zoo.krillDryWeight(:,:,wd,:) = repmat(zoo.krillDryWeight(:,:,1,:), [1 1 wn 1]) ./ wn;
zoo.krillWetWeight(:,:,wd,:) = repmat(zoo.krillWetWeight(:,:,1,:), [1 1 wn 1]) ./ wn;
zoo.krillC(:,:,wd,:) = repmat(zoo.krillC(:,:,1,:), [1 1 wn 1]) ./ wn;

% Organise the map_krill output. This should be a table saved as a .csv
% file that can be loaded into R.
lonmin = domain.longrid(1:end-1); lonmax = domain.longrid(2:end);
latmin = domain.latgrid(1:end-1); latmax = domain.latgrid(2:end);
lonmin = repmat(lonmin, [1 domain.nlat domain.nmonth]);
lonmax = repmat(lonmax, [1 domain.nlat domain.nmonth]);
latmin = repmat(reshape(latmin, 1, []), [domain.nlon 1 domain.nmonth]);
latmax = repmat(reshape(latmax, 1, []), [domain.nlon 1 domain.nmonth]);
month = repmat(reshape(domain.month, 1, 1, []), [domain.nlon domain.nlat 1]);
lonmin = lonmin(:); lonmax = lonmax(:);
latmin = latmin(:); latmax = latmax(:);
month = month(:);
value = map_krill(:);
map_krill_table = table(lonmin, lonmax, latmin, latmax, month, value);


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

[nlon, nlat, ~] = size(Dat.month);

% Filter the data by month
ind = ismember(Dat.month, months);
Dat = structfun(@(z) reshape(z(ind), nlon, nlat, []), Dat, ...
    'UniformOutput', false);
% Dat = structfun(@(z) z(ind), Dat, ...
%     'UniformOutput', false);

% Index data within mapped region
inmap = inpolygon(Dat.longitude, Dat.latitude, mv(1,:), mv(2,:));

% Omit data outside mapped region
% Dat = structfun(@(z) z(inmap), Dat, 'UniformOutput', false);

% Create struct for phytoplankton prey measurements -- similar to phys
prey = rmfield(phys, 'density');
prey.chl = zeros(domain.nlon, domain.nlat, domain.ndepth, domain.nmonth);

% Smooth the measurements into a regular grid to match map domain.
datgrid = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
for i = 1:domain.nlon
    indi = inmap & ...
        domain.longrid(i) < Dat.longitude & Dat.longitude <= domain.longrid(i+1);
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

map_chl = squeeze(datgrid); % output for plotting map of chlorophyll concentration

% Organise the map_chl output. This should be a table saved as a .csv
% file that can be loaded into R.
lonmin = domain.longrid(1:end-1); lonmax = domain.longrid(2:end);
latmin = domain.latgrid(1:end-1); latmax = domain.latgrid(2:end);
lonmin = repmat(lonmin, [1 domain.nlat domain.nmonth]);
lonmax = repmat(lonmax, [1 domain.nlat domain.nmonth]);
latmin = repmat(reshape(latmin, 1, []), [domain.nlon 1 domain.nmonth]);
latmax = repmat(reshape(latmax, 1, []), [domain.nlon 1 domain.nmonth]);
month = repmat(reshape(domain.month, 1, 1, []), [domain.nlon domain.nlat 1]);
lonmin = lonmin(:); lonmax = lonmax(:);
latmin = latmin(:); latmax = latmax(:);
month = month(:);
value = map_chl(:);
map_chl_table = table(lonmin, lonmax, latmin, latmax, month, value);

% Also output a high-resolution chlorophyll table that might be nice for
% the interactive map.
% The krill map uses a 9*3 degree lon*lat ratio for grid cells, so let's
% maintain that ratio for the chlorophyll.
switch makeHiResChlMap, case true
    dlon = domain.dlon / chlExtraResolution;
    dlat = domain.dlat / chlExtraResolution;
    longrid = domain.lon_range(1):dlon:domain.lon_range(2);
    latgrid = domain.lat_range(1):dlat:domain.lat_range(2);
    nlon = length(longrid) - 1; nlat = length(latgrid) - 1;
    % Find average values within each grid cell
    datgrid = nan(nlon, nlat, domain.nmonth);
    for i = 1:nlon
        indi = longrid(i) < Dat.longitude & Dat.longitude <= longrid(i+1);
        for j = 1:nlat
            indj = indi & ...
                latgrid(j) < Dat.latitude & Dat.latitude <= latgrid(j+1);
            for k = 1:domain.nmonth
                m = domain.month(k);
                indk = indj & Dat.month == m;
                if ~any(indk), continue; end
                datgrid(i,j,k) = mean(Dat.(vars_chl)(indk), 'omitnan');
            end
        end
%         disp(num2str(i / nlon))
    end
    lonmin = longrid(1:end-1); lonmax = longrid(2:end);
    latmin = latgrid(1:end-1); latmax = latgrid(2:end);
    lonmin = repmat(lonmin, [1 nlat domain.nmonth]);
    lonmax = repmat(lonmax, [1 nlat domain.nmonth]);
    latmin = repmat(reshape(latmin, 1, []), [nlon 1 domain.nmonth]);
    latmax = repmat(reshape(latmax, 1, []), [nlon 1 domain.nmonth]);
    month = repmat(reshape(domain.month, 1, 1, []), [nlon nlat 1]);
    lonmin = lonmin(:); lonmax = lonmax(:);
    latmin = latmin(:); latmax = latmax(:);
    month = month(:);
    value = datgrid(:);
    map_chl_table_highRes = table(lonmin, lonmax, latmin, latmax, month, value);
    switch saveOutput, case true
        path = fullfile(baseDirectory, 'MatLab', 'temp');
        if ~exist(path, 'dir')
            mkdir(fileparts(path), 'temp')
            addpath(genpath(path))
        end
        filename = ['chl_data_mapped_highRes_' num2str(round(dlon,2,'significant')) 'x' num2str(round(dlat,2,'significant')) '.csv'];
        filepath = fullfile(path, filename);
        writetable(map_chl_table_highRes, filepath)
    end
end











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
forc.krillWetWeight = zoo.krillWetWeight;
forc.krillDryWeight = zoo.krillDryWeight;
forc.krillC = zoo.krillC;
forc.noKrillData = noKrillData;
forc.krillLenMeanLog = zoo.krillLenMeanLog;
forc.krillLenSDLog = zoo.krillLenSDLog;
forc.krillDryWtMeanLog = zoo.krillDryWtMeanLog;
forc.krillDryWtSDLog = zoo.krillDryWtSDLog;
forc.krillWetWtMeanLog = zoo.krillWetWtMeanLog;
forc.krillWetWtSDLog = zoo.krillWetWtSDLog;
forc.krillCMeanLog = zoo.krillCMeanLog;
forc.krillCSDLog = zoo.krillCSDLog;

forcSize = [domain.nlon, domain.nlat, domain.ndepth, domain.nmonth];

correctSize = [all(size(forc.density_seawater) == forcSize), ...
    all(size(forc.chl_total) == forcSize), ... 
    all(size(forc.krill) == forcSize)];

if ~all(correctSize)
    warning('Some forcing data have incorrect dimension! Something is not right...')
end


%% Save
switch saveOutput, case true
    path = fullfile(baseDirectory, 'MatLab', 'temp');
    if ~exist(path, 'dir')
        mkdir(fileparts(path), 'temp')
        addpath(genpath(path))
    end
    filename = 'forcing.mat';
    filepath = fullfile(path, filename);
    save(filepath, 'forc')
    filename = 'domain.mat';
    filepath = fullfile(path, filename);
    save(filepath, 'domain')

    filename = 'krill_data_mapped.csv';
    filepath = fullfile(path, filename);
    writetable(map_krill_table, filepath)

    filename = 'chl_data_mapped.csv';
    filepath = fullfile(path, filename);
    writetable(map_chl_table, filepath)

    filename = 'chl_data_mapped_highRes.csv';
    filepath = fullfile(path, filename);
    writetable(map_chl_table_highRes, filepath)

end

%%
switch timeSetUp, case true
    setUpTime = toc / 60; fprintf('\n'); disp(append("finished: ", string(datetime('now')))); fprintf('\n');
    disp(['Forcing data set-up time: ' num2str(floor(setUpTime)) ' mins, ' ...
        num2str(floor(mod(60*setUpTime,60))) ' secs']); fprintf('\n\n')
end



