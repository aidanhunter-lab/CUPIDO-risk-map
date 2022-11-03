%% Load and organise forcing data
function [forc, domain] = prepare_forcing(baseDirectory, domain, varargin)

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
whichDepths = 1;
wd = false(domain.ndepth, 1); wd(whichDepths) = true;
zoo.krill(:,:,wd,:) = datgrid; % individuals/m^2
zoo.krill(:,:,~wd,:) = nan;

% Change units to number/m^3
zoo.krill =  zoo.krill .* ... 
    (0.5 * (domain.area(:,:,1:end-1) + domain.area(:,:,2:end))); % individuals/grid cell
zoo.krill =  zoo.krill ./ domain.volume; % individuals/m^3


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
    % length to weight
    try
        krillLength.WEIGHT_mg_dry = pars.W_dry_a .* krillLength.LENGTH_mm_ .^ pars.W_dry_b;
        krillLength.WEIGHT_mg_wet = pars.W_wet_a .* krillLength.LENGTH_mm_ .^ pars.W_wet_b;
        % carbon content from dry weight
        krillLength.CARBON_mg = pars.W_dry2c .* krillLength.WEIGHT_mg_dry;
    catch
        warning("Optional argument 'pars' is present but the allometric parameters are missing. See 'intialise_parameters.m'.")
    end
else
    filename = 'Unit Conversions.csv';
    tr = readtable(filename);
    % length to weight
    a = 'W_dry_a';
    b = 'W_dry_b';
    % tr.Unit(strcmp(tr.Parameter, a))
    % tr.Unit(strcmp(tr.Parameter, b))
    a = tr.Value(strcmp(tr.Parameter, a));
    b = tr.Value(strcmp(tr.Parameter, b));
    krillLength.WEIGHT_mg_dry = a .* krillLength.LENGTH_mm_ .^ b;
    a = 'W_wet_a';
    b = 'W_wet_b';
    a = tr.Value(strcmp(tr.Parameter, a));
    b = tr.Value(strcmp(tr.Parameter, b));
    krillLength.WEIGHT_mg_wet = a .* krillLength.LENGTH_mm_ .^ b;
    % carbon content from dry weight
    p = 'W_dry2c';
    s = 'Summer';
    p = tr.Value(strcmp(tr.Parameter, p) & strcmp(tr.Group, s));
    krillLength.CARBON_mg = p .* krillLength.WEIGHT_mg_dry;
end

% Match data to map grid - get weighted means and sd
% lengths
datgrid1 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
datgrid2 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
% dry weight
datgrid3 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
datgrid4 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
% wet weight
datgrid5 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
datgrid6 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
% carbon weight
datgrid7 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);
datgrid8 = nan(domain.nlon, domain.nlat, 1, domain.nmonth);

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

%             figure
%             h = histogram(dat.CARBON_mg(dat.EVENT == 15), 'FaceAlpha', 0.4)
%             hold on
%             histogram(dat.CARBON_mg(dat.EVENT == 17), 'FaceAlpha', 0.4)
%             histogram(dat.CARBON_mg(dat.EVENT == 18), 'FaceAlpha', 0.4)

            % Mean lengths weighted by event sample size
            ev = unique(dat.EVENT);
            evn = arrayfun(@(z) sum(dat.EVENT == z), ev); % number of measures per sample event
            % Summary stats on log scale
            meanLogLen = arrayfun(@(z) mean(log(dat.LENGTH_mm_(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogLen = sum(evn ./ sum(evn) .* meanLogLen);
            sdLogLen = arrayfun(@(z) std(log(dat.LENGTH_mm_(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogLen = sum(evn ./ sum(evn) .* sdLogLen);

            meanLogWdry = arrayfun(@(z) mean(log(dat.WEIGHT_mg_dry(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogWdry = sum(evn ./ sum(evn) .* meanLogWdry);
            sdLogWdry = arrayfun(@(z) std(log(dat.WEIGHT_mg_dry(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogWdry = sum(evn ./ sum(evn) .* sdLogWdry);

            meanLogWwet = arrayfun(@(z) mean(log(dat.WEIGHT_mg_wet(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogWwet = sum(evn ./ sum(evn) .* meanLogWwet);
            sdLogWwet = arrayfun(@(z) std(log(dat.WEIGHT_mg_wet(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogWwet = sum(evn ./ sum(evn) .* sdLogWwet);

            meanLogC = arrayfun(@(z) mean(log(dat.CARBON_mg(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            meanLogC = sum(evn ./ sum(evn) .* meanLogC);
            sdLogC = arrayfun(@(z) std(log(dat.CARBON_mg(dat.EVENT == z)), 'omitnan'), ev); % log-scale mean lengths
            sdLogC = sum(evn ./ sum(evn) .* sdLogC);

%             rr = normrnd(meanLogC, sdLogC, 1000, 1);
%             rr = exp(rr);
%             figure
%             histogram(rr)
            
            datgrid1(i,j,1,k) = meanLogLen;
            datgrid2(i,j,1,k) = sdLogLen;
            datgrid3(i,j,1,k) = meanLogWdry;
            datgrid4(i,j,1,k) = sdLogWdry;
            datgrid5(i,j,1,k) = meanLogWwet;
            datgrid6(i,j,1,k) = sdLogWwet;
            datgrid7(i,j,1,k) = meanLogC;
            datgrid8(i,j,1,k) = sdLogC;
        end
    end
end

A = nan(domain.nlon, domain.nlat, domain.ndepth, domain.nmonth);
zoo.krillLenMeanLog = A;
zoo.krillLenSDLog = A;
zoo.krillDryWtMeanLog = A;
zoo.krillDryWtSDLog = A;
zoo.krillWetWtMeanLog = A;
zoo.krillWetWtSDLog = A;
zoo.krillCMeanLog = A;
zoo.krillCSDLog = A;

% Distribute between modelled depth layers.
zoo.krillLenMeanLog(:,:,wd,:) = datgrid1;
zoo.krillLenSDLog(:,:,wd,:) = datgrid2;
zoo.krillDryWtMeanLog(:,:,wd,:) = datgrid3;
zoo.krillDryWtSDLog(:,:,wd,:) = datgrid4;
zoo.krillWetWtMeanLog(:,:,wd,:) = datgrid5;
zoo.krillWetWtSDLog(:,:,wd,:) = datgrid6;
zoo.krillCMeanLog(:,:,wd,:) = datgrid7;
zoo.krillCSDLog(:,:,wd,:) = datgrid8;

% Some density samples are not matched to length stats -- infill with
% nearest neighbours
gotL = ~isnan(zoo.krillLenMeanLog); % index data values
misL = ~isnan(zoo.krill) & ~gotL; % index data to infill
for j = 1:length(zoo.month)
    g = gotL(:,:,wd,j);
    i = misL(:,:,wd,j);
    li = find(i);
    [ri, ci] = find(i);
    fi = [ri, ci];
    [ri, ci] = find(g);
    gi = [ri, ci];
    ni = sum(i(:));
    for k = 1:ni
        % length
        ml = zoo.krillLenMeanLog(:,:,wd,j);
        ms = zoo.krillLenSDLog(:,:,wd,j);
        ii = fi(k,:);
        d = (sum(abs(gi - ii) .^ 2, 2)) .^ 0.5;
        md = find(d == min(d), 1);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillLenMeanLog(:,:,wd,j) = ml;
        zoo.krillLenSDLog(:,:,wd,j) = ms;
        % dry weight
        ml = zoo.krillDryWtMeanLog(:,:,wd,j);
        ms = zoo.krillDryWtSDLog(:,:,wd,j);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillDryWtMeanLog(:,:,wd,j) = ml;
        zoo.krillDryWtSDLog(:,:,wd,j) = ms;
        % wet weight
        ml = zoo.krillWetWtMeanLog(:,:,wd,j);
        ms = zoo.krillWetWtSDLog(:,:,wd,j);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillWetWtMeanLog(:,:,wd,j) = ml;
        zoo.krillWetWtSDLog(:,:,wd,j) = ms;
        % carbon weight
        ml = zoo.krillCMeanLog(:,:,wd,j);
        ms = zoo.krillCSDLog(:,:,wd,j);
        ml(li(k)) = ml(gi(md,1), gi(md,2));
        ms(li(k)) = ms(gi(md,1), gi(md,2));
        zoo.krillCMeanLog(:,:,wd,j) = ml;
        zoo.krillCSDLog(:,:,wd,j) = ms;
    end
end


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
end

%%
switch timeSetUp, case true
    setUpTime = toc / 60; fprintf('\n'); disp(append("finished: ", string(datetime('now')))); fprintf('\n');
    disp(['Forcing data set-up time: ' num2str(floor(setUpTime)) ' mins, ' ...
        num2str(floor(mod(60*setUpTime,60))) ' secs']); fprintf('\n\n')
end




