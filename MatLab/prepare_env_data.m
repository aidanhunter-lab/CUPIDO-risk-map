%% Prepare data for display on maps

%% Directories
thisFile = which('prepare_env_data.m');
dirBase = fileparts(fileparts(thisFile));
dirData = fullfile(dirBase, 'data');
dirMatLab = fullfile(dirBase, 'MatLab');
addpath(genpath(dirData))
addpath(genpath(dirMatLab))

%% Model domain
% Specify the modelled area and horizontal & vertical resolutions
domain = map_domain();
% Returns info on domain grid resolution [lon-lat coords] and grid cell
% areas & volumes [m^2 & m^3].

%% Prepare data
loadFromFile = false;
saveOutput = true;

%% Krill data

% Density (individuals / m2)
filename = 'krillbase_data.csv'; % Stored in ../data/KRILLBASE/density data/
krill = readtable(filename); % load data

% Years of krillbase data to retain: 2-element vector (min,max)
year_range_krillbase = [-inf,inf]; % Use all years
% Use day-time and/or night-time samples from krillbase
dayornight_krillbase = 'both'; % select 'day', 'night', or 'both'
if ~ismember(dayornight_krillbase, {'both','day','night'})
    error("Input argument dayornight_krillbase must be 'both', 'day', or 'night'.")
end
% Krillbase measurement variable to use as abundance value
var_krillbase = 'STANDARDISED_KRILL_UNDER_1M2'; % use the standardised measurement
% Depth layers with krill -- this is an artifact from previous modelling
% work where krill (numbers/m2) could be spread over multiple depth layers.
% Not relevant for this project => set to 1 (the surface layer).
krillDepths = 1;

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

% Filter by year
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

% Omit NaNs
krill = krill(~isnan(krill.(var_krillbase)),:);

% Create struct for krill measurements
zoo.lon = domain.lon;
zoo.lat = domain.lat;
zoo.depth = domain.depth;
zoo.month = domain.month;
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
%             v = mean(dat.(var_krillbase), 'omitnan');
            v = geomean(dat.(var_krillbase));
            datgrid(i,j,1,k) = v;
        end
    end
    percent_prog = floor(100 * i / domain.nlon);
    percent_prog_ = floor(100 * (i-1) / domain.nlon);
    if percent_prog ~= percent_prog_
        disp(['grid krill data ' num2str(percent_prog) '%'])
    end
end

% Index grid cells lacking data
noKrillData = isnan(datgrid);
% and set all these as NaN
% zoo.krill(repmat(noKrillData, [1, 1, domain.ndepth, 1])) = nan;

% Put all krill into top depth layer, later they may be distributed across depth
zoo.krill(:,:,1,:) = datgrid; % individuals/m^2
map_krill = squeeze(datgrid); % output for plotting map of krill abundance

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

switch saveOutput, case true
    filename = 'krill_data_mapped.csv';
    path = fullfile(dirData, 'KRILLBASE', 'compiled');
    if ~exist(path, 'dir')
        mkdir(path)
        addpath(genpath(path))
    end
    filepath = fullfile(path, filename);
    writetable(map_krill_table, filepath)
end


%% Phytoplankton data

filename = 'SeaWiFS_Phytoplankton-Size-Class-1997-2007-Southern-Ocean.mat'; % stored in ../data/chlorophyll/SeaWiFS/
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
% keepVars = [{'month','latitude','longitude'}, vars_chl];
keepVars = fields;
Dat = rmfield(Dat, fields(~ismember(fields, keepVars)));

[nlon, nlat, ~] = size(Dat.month);

% Filter the data by month
ind = ismember(Dat.month, domain.month);
Dat = structfun(@(z) reshape(z(ind), nlon, nlat, []), Dat, ...
    'UniformOutput', false);

% Index data within mapped region
inmap = inpolygon(Dat.longitude, Dat.latitude, mv(1,:), mv(2,:));

% Omit data outside mapped region
% Dat = structfun(@(z) z(inmap), Dat, 'UniformOutput', false);

% Create struct for phytoplankton prey measurements
prey.lon = domain.lon;
prey.lat = domain.lat;
prey.depth = domain.depth;
prey.month = domain.month;
A = nan(domain.nlon, domain.nlat, domain.ndepth, domain.nmonth);
prey.chl = A;

vars_chl = 'Tchl'; % select total chlorophyll from available data

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
            if ~any(indk(:)), continue; end
            dat = structfun(@(z) z(indk), Dat, 'UniformOutput', false);
%             v = mean(dat.(vars_chl), 'omitnan'); % this is not robust to multiple useVars values
            % The geometric mean is probably more appropriate
            vnan = isnan(dat.(vars_chl));
            vv = dat.(vars_chl)(~vnan);
            v = geomean(v);
            datgrid(i,j,1,k) = v;
        end
    end
    percent_prog = floor(100 * i / domain.nlon);
    percent_prog_ = floor(100 * (i-1) / domain.nlon);
    if percent_prog ~= percent_prog_
        disp(['grid 9x3 chl data ' num2str(percent_prog) '%'])
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

switch saveOutput, case true
    filename = ['chl_data_mapped_res_' num2str(round(domain.dlon,2,'significant')) 'x' num2str(round(domain.dlat,2,'significant')) '.csv'];
    path = fullfile(dirData, 'chlorophyll', 'SeaWiFS');
    if ~exist(path, 'dir')
        mkdir(path)
        addpath(genpath(path))
    end
    filepath = fullfile(path, filename);
    writetable(map_chl_table, filepath)
end

% Also output a high-resolution chlorophyll table that might be nice for
% the interactive map.
% The krill map uses a 9*3 degree lon*lat ratio for grid cells, so let's
% maintain that ratio for the chlorophyll.
chlExtraResolution = 3; % resolution increase factor

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
            if ~any(indk(:)), continue; end
            %                 datgrid(i,j,k) = mean(Dat.(vars_chl)(indk), 'omitnan');
            % Use geometric mean...
            v = Dat.(vars_chl)(indk);
            v = v(~isnan(v));
            datgrid(i,j,k) = geomean(v);
        end
    end
    percent_prog = floor(100 * i / nlon);
    percent_prog_ = floor(100 * (i-1) / nlon);
    if percent_prog ~= percent_prog_
        disp(['grid 3x1 chl data ' num2str(percent_prog) '%'])
    end
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
    filename = ['chl_data_mapped_res_' num2str(round(dlon,2,'significant')) 'x' num2str(round(dlat,2,'significant')) '.csv'];
    path = fullfile(dirData, 'chlorophyll', 'SeaWiFS');
    if ~exist(path, 'dir')
        mkdir(path)
        addpath(genpath(path))
    end
    filepath = fullfile(path, filename);
    writetable(map_chl_table_highRes, filepath)
end



%% SST

% The SST data are high resolution and span multiple decades, so the
% original data were roughly 0.5GB per month for all areas south of 45
% degrees. As the full data set would therefore require about 200 GB of
% storage, when downloading the data I took monthly and spatial means
% to reduce memory usage before saving the data. The data are saved in 2
% spatial resolutions: a 9*3 degree grid (matching the krill data) and a
% 3*1 degree grid (which works well with the interactive map). The 9*3 grid
% matches the spatial grid I intend to use for modelling.

filename_9x3 = 'sst_monthly_means_1982_2021_res_9x3.nc'; % stored in ../data/sst/ESA/
filename_3x1 = 'sst_monthly_means_1982_2021_res_3x1.nc';

ymin = 1982;
ymax = 2021;
% nyrs = ymax - ymin + 1;

% ncinfo(filename_3x1)
% ncdisp(filename_3x1)
% ncinfo(filename_9x3)
% ncdisp(filename_9x3)

lon_9x3 = ncread(filename_9x3, 'lon');
lat_9x3 = ncread(filename_9x3, 'lat');
time_9x3 = ncread(filename_9x3, 'time');
sst_9x3 = ncread(filename_9x3, 'analysed_sst');

lon_3x1 = ncread(filename_3x1, 'lon');
lat_3x1 = ncread(filename_3x1, 'lat');
time_3x1 = ncread(filename_3x1, 'time');
sst_3x1 = ncread(filename_3x1, 'analysed_sst');

nmonth_sst = length(time_9x3);

% Transform temperature from Kelvin to Celcius
K2C = @(x) x - 272.15;
sst_9x3 = K2C(sst_9x3);
sst_3x1 = K2C(sst_3x1);

% Convert data into table to save as .csv file.
% 1st 9x3, then 3x1 data
longrid = [lon_9x3 - 4.5; lon_9x3(end) + 4.5];
latgrid = [lat_9x3 - 1.5; lat_9x3(end) + 1.5];
% Trim excess data north of map domain boundary
latExcess = latgrid > max(domain.latgrid);
sst_9x3 = sst_9x3(:, ~latExcess(2:end), :);
latgrid = latgrid(latgrid <= max(domain.latgrid));

% Find min/max lat-lon values that define grid cells
lonmin = longrid(1:end-1); lonmax = longrid(2:end);
latmin = latgrid(1:end-1); latmax = latgrid(2:end);
lonmin = repmat(lonmin, [1 domain.nlat nmonth_sst]);
lonmax = repmat(lonmax, [1 domain.nlat nmonth_sst]);
latmin = repmat(reshape(latmin, 1, []), [domain.nlon 1 nmonth_sst]);
latmax = repmat(reshape(latmax, 1, []), [domain.nlon 1 nmonth_sst]);

year = ymin:ymax;
year = repmat(year, [12, 1]); year = year(:);
month = 1 + mod(time_9x3, 12);
month = repmat(reshape(month, 1, 1, []), [domain.nlon, domain.nlat, 1]);
year = repmat(reshape(year, 1, 1, []), [domain.nlon, domain.nlat, 1]);

lonmin = lonmin(:); lonmax = lonmax(:);
latmin = latmin(:); latmax = latmax(:);
month = month(:); year = year(:);
value = sst_9x3(:);
map_sst_table_9x3 = table(lonmin, lonmax, latmin, latmax, month, year, value);

% Now the 3x1 data
longrid = [lon_3x1 - 1.5; lon_3x1(end) + 1.5];
latgrid = [lat_3x1 - 0.5; lat_3x1(end) + 0.5];

% Trim excess data north of map domain boundary
latExcess = latgrid > max(domain.latgrid);
sst_3x1 = sst_3x1(:, ~latExcess(2:end), :);
latgrid = latgrid(~latExcess);
nlon_3x1 = length(longrid)-1;
nlat_3x1 = length(latgrid)-1;

% Find min/max lat-lon values that define grid cells
lonmin = longrid(1:end-1); lonmax = longrid(2:end);
latmin = latgrid(1:end-1); latmax = latgrid(2:end);
lonmin = repmat(lonmin, [1 nlat_3x1 nmonth_sst]);
lonmax = repmat(lonmax, [1 nlat_3x1 nmonth_sst]);
latmin = repmat(reshape(latmin, 1, []), [nlon_3x1 1 nmonth_sst]);
latmax = repmat(reshape(latmax, 1, []), [nlon_3x1 1 nmonth_sst]);

year = ymin:ymax;
year = repmat(year, [12, 1]); year = year(:);
month = 1 + mod(time_3x1, 12);
month = repmat(reshape(month, 1, 1, []), [nlon_3x1, nlat_3x1, 1]);
year = repmat(reshape(year, 1, 1, []), [nlon_3x1, nlat_3x1, 1]);

lonmin = lonmin(:); lonmax = lonmax(:);
latmin = latmin(:); latmax = latmax(:);
month = month(:); year = year(:);
value = sst_3x1(:);
map_sst_table_3x1 = table(lonmin, lonmax, latmin, latmax, month, year, value);


% Calculate SST anomoly...
% For each grid cell find the long-term monthly mean sst, spannning the
% whole time series.

% 1st the 9x3 data...
dat = map_sst_table_9x3; % coarse spatial resolution
yrs = sort(unique(dat.year));
nyrs = length(yrs);

a = table2array(dat);
d = unique(a(:,1:5), 'rows', 'stable');
n = size(d, 1); % number of unique combinations of grid cell and month
x = nan(n, nyrs); % temporary storage
for i = 1:nyrs
    y = yrs(i);
    t = a(:,6) == y;
    x(:,i) = a(t,7);
end
d(:,6) = mean(x, 2, 'omitnan'); % multi-year means for each grid cell and month
d = array2table(d, 'VariableNames', dat.Properties.VariableNames([1:5 7]));

% d_9x3 = d;

% Find SST anomoly as difference between measurements and long term mean

d.group = (1:size(d, 1))';
d_ = d; d_.value = [];
dat = join(dat, d_);

x = nan(max(dat.group), nyrs);

for i = 1:nyrs
    y = yrs(i);
    t = dat.year == y;
    dy = dat(t,:);
    v = dy.value - d.value;
    x(:,i) = v;
end

map_sst_table_9x3.anomoly = x(:);

% Now repeat for the finer res 3x1 data

dat = map_sst_table_3x1;
yrs = sort(unique(dat.year));
nyrs = length(yrs);

a = table2array(dat);
d = unique(a(:,1:5), 'rows', 'stable');
n = size(d, 1); % number of unique combinations of grid cell and month
x = nan(n, nyrs); % temporary storage
for i = 1:nyrs
    y = yrs(i);
    t = a(:,6) == y;
    x(:,i) = a(t,7);
end
d(:,6) = mean(x, 2, 'omitnan'); % multi-year means for each grid cell and month
d = array2table(d, 'VariableNames', dat.Properties.VariableNames([1:5 7]));

% d_3x1 = d;

% Find SST anomoly as difference between measurements and long term mean

d.group = (1:size(d, 1))';
d_ = d; d_.value = [];
dat = join(dat, d_);

x = nan(max(dat.group), nyrs);

for i = 1:nyrs
    y = yrs(i);
    t = dat.year == y;
    dy = dat(t,:);
    v = dy.value - d.value;
    x(:,i) = v;
end

map_sst_table_3x1.anomoly = x(:);


% Calculate SST long-term linear trend
% For each grid cell and month fit a linear model of SST vs year spanning
% the whole time series. Do the same for aggregation of summer months (Jan-Mar)

% 1st the 9x3 data...
dat = map_sst_table_9x3;

d = unique(dat(:,1:4), 'stable');
d.grid_cell = (1:height(d))';
dat = join(dat, d);
ng = max(dat.grid_cell);
hd = height(dat);
Months = datestr(datetime(1,1:12,1),'mmm');
dat.linearTrend_Summer = nan(hd, 1);
dat.pValue_Summer = nan(hd, 1);
for i = 1:12
    dat.(['linearTrend_' Months(i,:)]) = nan(hd, 1);
    dat.(['pValue_' Months(i,:)]) = nan(hd, 1);
end
for i = 1:ng
    ii = dat.grid_cell == i;
    y = dat(ii,:);
    if all(isnan(y.value)), continue; end
    a = y(ismember(y.month, 1:3),:); % extract summer months
    av = reshape(a.value, 3, []);
    av = mean(av, 'omitnan'); av = av(:); % average over summer months    
    ma = fitlm(yrs, av); % fit a linear model
    tra = ma.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pva = ma.Coefficients.pValue(2); % get p-values to gauge significance
    dat.linearTrend_Summer(ii) = tra;
    dat.pValue_Summer(ii) = pva;
    % repeat for each month independently
    for j = 1:12
        Mon = Months(j,:);
        yj = y(y.month == j,:);
        m = fitlm(yj.year, yj.value);
        tr = m.Coefficients.Estimate(2);
        pv = m.Coefficients.pValue(2);
        dat.(['linearTrend_' Mon])(ii) = tr;
        dat.(['pValue_' Mon])(ii) = pv;
    end
    percent_prog = floor(100 * i / ng);
    percent_prog_ = floor(100 * (i-1) / ng);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 9x3 sst data ' num2str(percent_prog) '%'])
    end
end

dat.grid_cell = [];

map_sst_table_9x3 = dat;

% repeat for the higher res data
dat = map_sst_table_3x1;

d = unique(dat(:,1:4), 'stable');
d.grid_cell = (1:height(d))';
dat = join(dat, d);
ng = max(dat.grid_cell);
hd = height(dat);
dat.linearTrend_Summer = nan(hd, 1);
dat.pValue_Summer = nan(hd, 1);
for i = 1:12
    dat.(['linearTrend_' Months(i,:)]) = nan(hd, 1);
    dat.(['pValue_' Months(i,:)]) = nan(hd, 1);
end
for i = 1:ng
    disp([num2str(round(100 * i / ng, 2)) '%'])
    ii = dat.grid_cell == i;
    y = dat(ii,:);
    if all(isnan(y.value)), continue; end
    a = y(ismember(y.month, 1:3),:); % extract summer months
    av = reshape(a.value, 3, []);
    av = mean(av, 'omitnan'); av = av(:); % average over summer months    
    ma = fitlm(yrs, av); % fit a linear model
    tra = ma.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pva = ma.Coefficients.pValue(2); % get p-values to gauge significance
    dat.linearTrend_Summer(ii) = tra;
    dat.pValue_Summer(ii) = pva;
    % repeat for each month independently
    for j = 1:12
        Mon = Months(j,:);
        yj = y(y.month == j,:);
        m = fitlm(yj.year, yj.value);
        tr = m.Coefficients.Estimate(2);
        pv = m.Coefficients.pValue(2);
        dat.(['linearTrend_' Mon])(ii) = tr;
        dat.(['pValue_' Mon])(ii) = pv;
    end
    percent_prog = floor(100 * i / ng);
    percent_prog_ = floor(100 * (i-1) / ng);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 3x1 sst data ' num2str(percent_prog) '%'])
    end

end

dat.grid_cell = [];
map_sst_table_3x1 = dat;

% Filter this data before saving. I need to reduce memory requiremenrs for
% a speedy interactive map...

% For the trend data, all we need are single values for each grid cell --
% also include the p-values
coords_9x3 = unique(map_sst_table_9x3(:,1:4), 'stable');
sst_trend_9x3 = map_sst_table_9x3;
sst_trend_9x3 = sst_trend_9x3(:,[1:4 9:end]);
sst_trend_9x3 = sst_trend_9x3(~all(isnan(table2array(sst_trend_9x3(:,5:end))), 2),:);
sst_trend_9x3 = unique(sst_trend_9x3, 'stable');
sst_trend_9x3 = outerjoin(coords_9x3, sst_trend_9x3);
sst_trend_9x3 = sst_trend_9x3(:,[1:4 9:end]);
sst_trend_9x3.Properties.VariableNames = strrep(sst_trend_9x3.Properties.VariableNames, '_coords_9x3', '');
sst_trend_9x3.Properties.VariableNames = strrep(sst_trend_9x3.Properties.VariableNames, 'linearTrend_', '');

coords_3x1 = unique(map_sst_table_3x1(:,1:4), 'stable');
sst_trend_3x1 = map_sst_table_3x1;
sst_trend_3x1 = sst_trend_3x1(:,[1:4 9:end]);
sst_trend_3x1 = sst_trend_3x1(~all(isnan(table2array(sst_trend_3x1(:,5:end))), 2),:);
sst_trend_3x1 = unique(sst_trend_3x1, 'stable');
sst_trend_3x1 = outerjoin(coords_3x1, sst_trend_3x1);
sst_trend_3x1 = sst_trend_3x1(:,[1:4 9:end]);
sst_trend_3x1.Properties.VariableNames = strrep(sst_trend_3x1.Properties.VariableNames, '_coords_3x1', '');
sst_trend_3x1.Properties.VariableNames = strrep(sst_trend_3x1.Properties.VariableNames, 'linearTrend_', '');

% The temperature anomaly varies month to month across all years and grid
% cells
sst_anomaly_9x3 = map_sst_table_9x3(:,1:8);
sst_anomaly_3x1 = map_sst_table_3x1(:,1:8);

switch saveOutput, case true
    path = fullfile(dirData, 'sst', 'ESA', 'compiled');
    if ~exist(path, 'dir')
        mkdir(path)
        addpath(genpath(path))
    end

    filename = 'sst_trend_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(sst_trend_9x3, filepath)

    filename = 'sst_trend_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(sst_trend_3x1, filepath)

    filename = 'sst_anomaly_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(sst_anomaly_9x3, filepath)

    filename = 'sst_anomaly_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(sst_anomaly_3x1, filepath)
end



% Try the 5-step time series method described in the Global Ocean Sea Surface
% Temperature trend map from Observations Reprocessing.
% https://data.marine.copernicus.eu/product/GLOBAL_OMI_TEMPSAL_sst_trend/description

% 1. The daily analyses were averaged to create monthly means.
% 2. A climatology was calculated by averaging the monthly means over the 
% period 1993 - 2014.
% 3. Monthly anomalies were calculated by differencing the monthly means and
% the climatology.
% 4. The time series for each grid cell was passed through the X11 seasonal
% adjustment procedure, which decomposes a time series into a residual
% seasonal component, a trend component and errors (e.g., Pezzulli et al., 2005).
% The trend component is a filtered version of the monthly time series.
% 5. The slope of the trend component was calculated using a robust method
% (Sen 1968). The method also calculates the 95% confidence range in the slope.

% The monthly mean values (1) are already tabled.
% Steps 2 & 3, the climatologies and monthly anomalies, is similar to the
% anomalies calculated above. Repeat this, using the 1993-2014 period, to 
% harmonise my outputs with the existing trend map.
% Step 4 is the main difference -- time series method to account for
% seasonality.
% Step 5 is similar to the linear modelling above -- estimating the
% regression coefficient.

% Work with the coarsely resolved data
dat = map_sst_table_9x3(:,1:7);
% Index grid cell by a single column
d = unique(dat(:,1:4));
ngridcells = height(d);
d.gridCell = (1:ngridcells)';
dat = join(dat,d);
% Sort data
[~, I] = sortrows([dat.gridCell, dat.year, dat.month]);
dat = dat(I,:);
% Calculate the climatologies as monthly means over 1993-2014
climatology = nan(12, ngridcells);
climWindow = [1993, 2014];
ind = climWindow(1) <= dat.year & dat.year <= climWindow(2);
for i = 1:ngridcells
    d = dat(ind & dat.gridCell == i,:);
    for j = 1:12
        dm = d(d.month == j,:);
        climatology(j,i) = mean(dm.value, 'omitnan');
    end
end
% Calculate anomalies as difference from climatology
climatology3d = repmat(reshape(climatology, 12, 1, ngridcells), [1, nyrs, 1]);
dat.anomaly = dat.value - climatology3d(:);
% Time series analysis using X11 algorithm
dat.seasonal_factor = nan(height(dat), 1);
dat.random_factor = nan(height(dat), 1);
dat.trend = nan(height(dat), 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    vx = x11(d.anomaly, 12, 'add');
    dat.seasonal_factor(ind) = vx.d10;
    dat.random_factor(ind) = vx.d13;
    dat.trend(ind) = vx.d12;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
            disp(['decompose 9x3 sst data ' num2str(percent_prog) '%'])
    end
end

% Plot time series decomposition outputs
plot_time_series_decomposition = false;
switch plot_time_series_decomposition, case true
    gc = 520; % choose grid cell
    ind = dat.gridCell == gc;
    subplot(2,2,1)
    plot(dat.anomaly(ind))
    title('temperature anomaly')
    subplot(2,2,2)
    plot(dat.seasonal_factor(ind))
    title('seasonal factor')
    subplot(2,2,3)
    plot(dat.random_factor(ind))
    title('random factor')
    subplot(2,2,4)
    plot(dat.trend(ind))
    title('trend')
    sgtitle(['grid cell ' num2str(gc)])
end

% Calculate the linear slope of the trend.
% This will produce a single value per grid cell, unlike my previous method
% of calculating separate values for each month and grid cell.
rateOfChange = nan(ngridcells, 1);
pvals = nan(ngridcells, 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    y = d.trend; % temperature anomaly corrected for seasonality
    if all(isnan(y)), continue; end
    x = (1:length(y))' / 12; % independent variable: time (years)
    model = fitlm(x, y); % fit a linear model
    coef = model.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pval = model.Coefficients.pValue(2); % get p-values to gauge significance
    rateOfChange(i) = coef;
    pvals(i) = pval;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 9x3 sst data ' num2str(percent_prog) '%'])
    end
end

d = unique(dat(:,{'lonmin' 'lonmax', 'latmin' 'latmax' 'gridCell'}));
[~,I] = sort(d.gridCell);
d = d(I,:);
d.linearTrend = rateOfChange;
d.pValue = pvals;
dat_trend_9x3 = d;

% Now repeat for higer res data
dat = map_sst_table_3x1(:,1:7);
% Index grid cell by a single column
d = unique(dat(:,1:4));
ngridcells = height(d);
d.gridCell = (1:ngridcells)';
dat = join(dat,d);
% Sort data
[~, I] = sortrows([dat.gridCell, dat.year, dat.month]);
dat = dat(I,:);
% Calculate the climatologies as monthly means over 1993-2014
climatology = nan(12, ngridcells);
climWindow = [1993, 2014];
ind = climWindow(1) <= dat.year & dat.year <= climWindow(2);
for i = 1:ngridcells
    d = dat(ind & dat.gridCell == i,:);
    for j = 1:12
        dm = d(d.month == j,:);
        climatology(j,i) = mean(dm.value, 'omitnan');
    end
end
% Calculate anomalies as difference from climatology
climatology3d = repmat(reshape(climatology, 12, 1, ngridcells), [1, nyrs, 1]);
dat.anomaly = dat.value - climatology3d(:);
% Time series analysis using X11 algorithm
dat.seasonal_factor = nan(height(dat), 1);
dat.random_factor = nan(height(dat), 1);
dat.trend = nan(height(dat), 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    vx = x11(d.anomaly, 12, 'add');
    dat.seasonal_factor(ind) = vx.d10;
    dat.random_factor(ind) = vx.d13;
    dat.trend(ind) = vx.d12;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['decompose 3x1 sst data ' num2str(percent_prog) '%'])
    end
end
% Plot time series decomposition outputs
switch plot_time_series_decomposition, case true
    gc = 1000; % choose grid cell
    ind = dat.gridCell == gc;
    subplot(2,2,1)
    plot(dat.anomaly(ind))
    title('temperature anomaly')
    subplot(2,2,2)
    plot(dat.seasonal_factor(ind))
    title('seasonal factor')
    subplot(2,2,3)
    plot(dat.random_factor(ind))
    title('random factor')
    subplot(2,2,4)
    plot(dat.trend(ind))
    title('trend')
    sgtitle(['grid cell ' num2str(gc)])
end
% Calculate the linear slope of the trend.
% This will produce a single value per grid cell, unlike my previous method
% of calculating separate values for each month and grid cell.
rateOfChange = nan(ngridcells, 1);
    pvals = nan(ngridcells, 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    y = d.trend; % temperature anomaly corrected for seasonality
    if all(isnan(y)), continue; end
    x = (1:length(y))' / 12; % independent variable: time (years)
    model = fitlm(x, y); % fit a linear model
    coef = model.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pval = model.Coefficients.pValue(2); % get p-values to gauge significance
    rateOfChange(i) = coef;
    pvals(i) = pval;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 3x1 sst data ' num2str(percent_prog) '%'])
    end
end

d = unique(dat(:,{'lonmin' 'lonmax', 'latmin' 'latmax' 'gridCell'}));
[~,I] = sort(d.gridCell);
d = d(I,:);
d.linearTrend = rateOfChange;
d.pValue = pvals;
dat_trend_3x1 = d;


dat_trend_9x3.gridCell = [];
dat_trend_3x1.gridCell = [];

switch saveOutput, case true
    path = fullfile(dirData, 'sst', 'ESA', 'compiled');
    if ~exist(path, 'dir')
        mkdir(path)
        addpath(genpath(path))
    end

    filename = 'sst_time_series_trend_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(dat_trend_9x3, filepath)

    filename = 'sst_time_series_trend_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(dat_trend_3x1, filepath)
end


%% pH

% The pH data are monthly estimates with a 1x1 degree lat-lon resolution.

filename = 'pH_19850115_20211215.nc'; % stored in ../data/pH/SOCAT/

ncinfo(filename)
ncdisp(filename)

lon = ncread(filename, 'longitude');
lat = ncread(filename, 'latitude');
time = ncread(filename, 'time');
pH = ncread(filename, 'ph');
pH_sd = ncread(filename, 'ph_uncertainty');

% pH_orig = pH;
% pH_sd_orig = pH_sd;

% Convert the time variable to MatLab standard time
pH_startDate = [1950 1 1 0 0 0];
adjustTime = datenum(pH_startDate);
time = time / 24; % days since 1/1/1950
time = time + adjustTime;
ntime = length(time);

% Convert lat-lons to degrees N (-90,90) and degrees E (-180,180)
lon(lon > 180) = lon(lon > 180) - 360;
[~,I] = sort(lon);
lon = lon(I);
pH = pH(I,:,:);
pH_sd = pH_sd(I,:,:);

% Omit data outside mapped region
ilon = domain.lon_range(1) <= lon & lon <= domain.lon_range(2);
ilat = domain.lat_range(1) <= lat & lat <= domain.lat_range(2);
lon = lon(ilon);
lat = lat(ilat);
% nlon = length(lon);
% nlat = length(lat);
pH = pH(ilon,ilat,:);
pH_sd = pH_sd(ilon,ilat,:);

% Take averages to smooth data into required spatial resolution(s).
datgrid = nan(domain.nlon, domain.nlat, ntime);
pH_grid = datgrid;
pH_sd_grid = datgrid;
for i = 1:domain.nlon
    indi =  domain.longrid(i) <= lon & lon <= domain.longrid(i+1);
    if ~any(indi), continue; end
    for j = 1:domain.nlat
        indj = domain.latgrid(j) <= lat & lat <= domain.latgrid(j+1);
        if ~any(indj), continue; end
        d = pH(indi, indj, :);
        d = mean(d, [1, 2], 'omitnan');
        pH_grid(i,j,:) = d(:);
        d = pH_sd(indi, indj, :);
        d = mean(d .^ 2, [1, 2], 'omitnan') .^ 0.5;
        pH_sd_grid(i,j,:) = d(:);
    end
    percent_prog = floor(100 * i / domain.nlon);
    percent_prog_ = floor(100 * (i-1) / domain.nlon);
    if percent_prog ~= percent_prog_
        disp(['grid 9x3 pH data ' num2str(percent_prog) '%'])
    end
end

% Convert data into table to save as .csv file.
longrid = domain.longrid;
latgrid = domain.latgrid;
lonmin = longrid(1:end-1); lonmax = longrid(2:end);
latmin = latgrid(1:end-1); latmax = latgrid(2:end);
lonmin = repmat(lonmin, [1 domain.nlat ntime]);
lonmax = repmat(lonmax, [1 domain.nlat ntime]);
latmin = repmat(reshape(latmin, 1, []), [domain.nlon 1 ntime]);
latmax = repmat(reshape(latmax, 1, []), [domain.nlon 1 ntime]);

[year, month] = datevec(double(time));
year = repmat(reshape(year, 1, 1, []), [domain.nlon, domain.nlat, 1]);
month = repmat(reshape(month, 1, 1, []), [domain.nlon, domain.nlat, 1]);

lonmin = lonmin(:); lonmax = lonmax(:);
latmin = latmin(:); latmax = latmax(:);
year = year(:); month = month(:);
pH = pH_grid(:);
pH_sd = pH_sd_grid(:);

map_pH_table_9x3 = table(lonmin, lonmax, latmin, latmax, month, year, pH, pH_sd);

% Create a higher resolution (3x1) data table for interactive map

lon = ncread(filename, 'longitude');
lat = ncread(filename, 'latitude');
time = ncread(filename, 'time');
pH = ncread(filename, 'ph');
pH_sd = ncread(filename, 'ph_uncertainty');

% Convert the time variable to MatLab standard time
pH_startDate = [1950 1 1 0 0 0];
adjustTime = datenum(pH_startDate);
time = time / 24; % days since 1/1/1950
time = time + adjustTime;
ntime = length(time);

% Convert lat-lons to degrees N (-90,90) and degrees E (-180,180)
lon(lon > 180) = lon(lon > 180) - 360;
[~,I] = sort(lon);
lon = lon(I);
pH = pH(I,:,:);
pH_sd = pH_sd(I,:,:);

% Omit data outside mapped region
ilon = domain.lon_range(1) <= lon & lon <= domain.lon_range(2);
ilat = domain.lat_range(1) <= lat & lat <= domain.lat_range(2);
lon = lon(ilon);
lat = lat(ilat);
% nlon = length(lon);
% nlat = length(lat);
pH = pH(ilon,ilat,:);
pH_sd = pH_sd(ilon,ilat,:);

% Take averages to smooth data into required spatial resolution(s).
longrid = domain.lon_range(1):3:domain.lon_range(2);
longrid = longrid(:);
latgrid = domain.lat_range(1):1:domain.lat_range(2);
nlon = length(longrid) - 1;
nlat = length(latgrid) - 1;
datgrid = nan(nlon, nlat, ntime);
pH_grid = datgrid;
pH_sd_grid = datgrid;
for i = 1:nlon
    indi =  longrid(i) <= lon & lon <= longrid(i+1);
    if ~any(indi), continue; end
    for j = 1:nlat
        indj = latgrid(j) <= lat & lat <= latgrid(j+1);
        if ~any(indj), continue; end
        d = pH(indi, indj, :);
        d = mean(d, [1, 2], 'omitnan');
        pH_grid(i,j,:) = d(:);
        d = pH_sd(indi, indj, :);
        d = mean(d .^ 2, [1, 2], 'omitnan') .^ 0.5;
        pH_sd_grid(i,j,:) = d(:);
    end
    percent_prog = floor(100 * i / nlon);
    percent_prog_ = floor(100 * (i-1) / nlon);
    if percent_prog ~= percent_prog_
        disp(['grid 3x1 pH data ' num2str(percent_prog) '%'])
    end
end

% Convert data into table to save as .csv file.
lonmin = longrid(1:end-1); lonmax = longrid(2:end);
latmin = latgrid(1:end-1); latmax = latgrid(2:end);
lonmin = repmat(lonmin, [1 nlat ntime]);
lonmax = repmat(lonmax, [1 nlat ntime]);
latmin = repmat(latmin, [nlon 1 ntime]);
latmax = repmat(latmax, [nlon 1 ntime]);

[year, month] = datevec(double(time));
year = repmat(reshape(year, 1, 1, []), [nlon, nlat, 1]);
month = repmat(reshape(month, 1, 1, []), [nlon, nlat, 1]);

lonmin = lonmin(:); lonmax = lonmax(:);
latmin = latmin(:); latmax = latmax(:);
year = year(:); month = month(:);
pH = pH_grid(:);
pH_sd = pH_sd_grid(:);

map_pH_table_3x1 = table(lonmin, lonmax, latmin, latmax, month, year, pH, pH_sd);

% Calculate pH anomoly...
% For each grid cell find the long-term monthly mean pH, spannning the
% whole time series.

% 1st the 9x3 data...
dat = map_pH_table_9x3;
yrs = sort(unique(dat.year));
nyrs = length(yrs);

a = table2array(dat);
d = unique(a(:,1:5), 'rows', 'stable');
n = size(d, 1); % number of unique combinations of grid cell and month
x = nan(n, nyrs); % temporary storage
for i = 1:nyrs
    y = yrs(i);
    t = a(:,6) == y;
    x(:,i) = a(t,7);
end
d(:,6) = mean(x, 2, 'omitnan'); % multi-year means for each grid cell and month
d = array2table(d, 'VariableNames', dat.Properties.VariableNames([1:5 7]));

d_9x3 = d;

% Find pH anomoly as difference between measurements and long term mean

d.group = (1:size(d, 1))';
d_ = d; d_.pH = [];
dat = join(dat, d_);

x = nan(max(dat.group), nyrs);

for i = 1:nyrs
    y = yrs(i);
    t = dat.year == y;
    dy = dat(t,:);
    v = dy.pH - d.pH;
    x(:,i) = v;
end

map_pH_table_9x3.anomaly = x(:);

% Now repeat for the finer res 3x1 data

dat = map_pH_table_3x1;
yrs = sort(unique(dat.year));
nyrs = length(yrs);

a = table2array(dat);
d = unique(a(:,1:5), 'rows', 'stable');
n = size(d, 1); % number of unique combinations of grid cell and month
x = nan(n, nyrs); % temporary storage
for i = 1:nyrs
    y = yrs(i);
    t = a(:,6) == y;
    x(:,i) = a(t,7);
end
d(:,6) = mean(x, 2, 'omitnan'); % multi-year means for each grid cell and month
d = array2table(d, 'VariableNames', dat.Properties.VariableNames([1:5 7]));

d_3x1 = d;

% Find pH anomoly as difference between measurements and long term mean

d.group = (1:size(d, 1))';
d_ = d; d_.pH = [];
dat = join(dat, d_);

x = nan(max(dat.group), nyrs);

for i = 1:nyrs
    y = yrs(i);
    t = dat.year == y;
    dy = dat(t,:);
    v = dy.pH - d.pH;
    x(:,i) = v;
end

map_pH_table_3x1.anomaly = x(:);

% Calculate pH long-term linear trend
% For each grid cell and month fit a linear model of pH vs year spanning
% the whole time series. Do the same for aggregation of summer months (Jan-Mar)

% 1st the 9x3 data...
dat = map_pH_table_9x3;

d = unique(dat(:,1:4), 'stable');
d.grid_cell = (1:height(d))';
dat = join(dat, d);
ng = max(dat.grid_cell);
hd = height(dat);
Months = datestr(datetime(1,1:12,1),'mmm');
dat.linearTrend_Summer = nan(hd, 1);
dat.pValue_Summer = nan(hd, 1);
for i = 1:12
    dat.(['linearTrend_' Months(i,:)]) = nan(hd, 1);
    dat.(['pValue_' Months(i,:)]) = nan(hd, 1);
end
for i = 1:ng
    ii = dat.grid_cell == i;
    y = dat(ii,:);
    if all(isnan(y.pH)), continue; end
    a = y(ismember(y.month, 1:3),:); % extract summer months
    av = reshape(a.pH, 3, []);
    av = mean(av, 'omitnan'); av = av(:); % average over summer months    
    ma = fitlm(yrs, av); % fit a linear model
    tra = ma.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pva = ma.Coefficients.pValue(2); % get p-values to gauge significance
    dat.linearTrend_Summer(ii) = tra;
    dat.pValue_Summer(ii) = pva;
    % repeat for each month independently
    for j = 1:12
%         disp(num2str(j))
        Mon = Months(j,:);
        yj = y(y.month == j,:);
        if sum(~isnan(yj.pH)) < 2, continue; end % skip if too few values for linear model (I suspect interannual variability in ice cover may be an issue here)
        m = fitlm(yj.year, yj.pH);
        tr = m.Coefficients.Estimate(2);
        pv = m.Coefficients.pValue(2);
        dat.(['linearTrend_' Mon])(ii) = tr;
        dat.(['pValue_' Mon])(ii) = pv;
    end
    percent_prog = floor(100 * i / ng);
    percent_prog_ = floor(100 * (i-1) / ng);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 9x3 pH data ' num2str(percent_prog) '%'])
    end
end

dat.grid_cell = [];

map_pH_table_9x3 = dat;

% repeat for the higher res data

dat = map_pH_table_3x1;

d = unique(dat(:,1:4), 'stable');
d.grid_cell = (1:height(d))';
dat = join(dat, d);
ng = max(dat.grid_cell);
hd = height(dat);
Months = datestr(datetime(1,1:12,1),'mmm');
dat.linearTrend_Summer = nan(hd, 1);
dat.pValue_Summer = nan(hd, 1);
for i = 1:12
    dat.(['linearTrend_' Months(i,:)]) = nan(hd, 1);
    dat.(['pValue_' Months(i,:)]) = nan(hd, 1);
end
for i = 1:ng
    ii = dat.grid_cell == i;
    y = dat(ii,:);
    if all(isnan(y.pH)), continue; end
    a = y(ismember(y.month, 1:3),:); % extract summer months
    av = reshape(a.pH, 3, []);
    av = mean(av, 'omitnan'); av = av(:); % average over summer months    
    ma = fitlm(yrs, av); % fit a linear model
    tra = ma.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pva = ma.Coefficients.pValue(2); % get p-values to gauge significance
    dat.linearTrend_Summer(ii) = tra;
    dat.pValue_Summer(ii) = pva;
    % repeat for each month independently
    for j = 1:12
%         disp(num2str(j))
        Mon = Months(j,:);
        yj = y(y.month == j,:);
        if sum(~isnan(yj.pH)) < 2, continue; end % skip if too few values for linear model (I suspect interannual variability in ice cover may be an issue here)
        m = fitlm(yj.year, yj.pH);
        tr = m.Coefficients.Estimate(2);
        pv = m.Coefficients.pValue(2);
        dat.(['linearTrend_' Mon])(ii) = tr;
        dat.(['pValue_' Mon])(ii) = pv;
    end
    percent_prog = floor(100 * i / ng);
    percent_prog_ = floor(100 * (i-1) / ng);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 3x1 pH data ' num2str(percent_prog) '%'])
    end
end

dat.grid_cell = [];

map_pH_table_3x1 = dat;


% Filter this data before saving. Reduce memory requiremenrs for a speedy
% interactive map...

% For the trend data, all we need are single values for each grid cell
coords_9x3 = unique(map_pH_table_9x3(:,1:4), 'stable');
pH_trend_9x3 = map_pH_table_9x3;
pH_trend_9x3 = pH_trend_9x3(:,[1:4 10:end]);
pH_trend_9x3 = pH_trend_9x3(~all(isnan(table2array(pH_trend_9x3(:,5:end))), 2),:);
a = table2array(pH_trend_9x3);
a(isnan(a)) = inf; % replace NaNs so that unique works 'properly'
a = unique(a, 'rows', 'stable');
a(a == inf) = nan;
pH_trend_9x3 = array2table(a, 'VariableNames', pH_trend_9x3.Properties.VariableNames);
pH_trend_9x3 = outerjoin(coords_9x3, pH_trend_9x3);
pH_trend_9x3 = pH_trend_9x3(:,[1:4 9:end]);
pH_trend_9x3.Properties.VariableNames = strrep(pH_trend_9x3.Properties.VariableNames, '_coords_9x3', '');
pH_trend_9x3.Properties.VariableNames = strrep(pH_trend_9x3.Properties.VariableNames, 'linearTrend_', '');

coords_3x1 = unique(map_pH_table_3x1(:,1:4), 'stable');
pH_trend_3x1 = map_pH_table_3x1;
pH_trend_3x1 = pH_trend_3x1(:,[1:4 10:end]);
pH_trend_3x1 = pH_trend_3x1(~all(isnan(table2array(pH_trend_3x1(:,5:end))), 2),:);
a = table2array(pH_trend_3x1);
a(isnan(a)) = inf; % replace NaNs so that unique works 'properly'
a = unique(a, 'rows', 'stable');
a(a == inf) = nan;
pH_trend_3x1 = array2table(a, 'VariableNames', pH_trend_3x1.Properties.VariableNames);
pH_trend_3x1 = outerjoin(coords_3x1, pH_trend_3x1);
pH_trend_3x1 = pH_trend_3x1(:,[1:4 9:end]);
pH_trend_3x1.Properties.VariableNames = strrep(pH_trend_3x1.Properties.VariableNames, '_coords_3x1', '');
pH_trend_3x1.Properties.VariableNames = strrep(pH_trend_3x1.Properties.VariableNames, 'linearTrend_', '');

% The pH anomaly varies month to month across all years and grid cells
pH_anomaly_9x3 = map_pH_table_9x3(:,1:9);
pH_anomaly_3x1 = map_pH_table_3x1(:,1:9);


switch saveOutput, case true
    path = fullfile(dirData, 'pH', 'SOCAT', 'compiled');
    if ~exist(path, 'dir')
        mkdir(path)
        addpath(genpath(path))
    end

    filename = 'pH_trend_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(pH_trend_9x3, filepath)

    filename = 'pH_trend_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(pH_trend_3x1, filepath)

    filename = 'pH_anomaly_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(pH_anomaly_9x3, filepath)

    filename = 'pH_anomaly_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(pH_anomaly_3x1, filepath)

end



% Try the 5-step time series method described in the Global Ocean Sea Surface
% Temperature trend map from Observations Reprocessing.
% https://data.marine.copernicus.eu/product/GLOBAL_OMI_TEMPSAL_sst_trend/description

% 1. The daily analyses were averaged to create monthly means.
% 2. A climatology was calculated by averaging the monthly means over the 
% period 1993 - 2014.
% 3. Monthly anomalies were calculated by differencing the monthly means and
% the climatology.
% 4. The time series for each grid cell was passed through the X11 seasonal
% adjustment procedure, which decomposes a time series into a residual
% seasonal component, a trend component and errors (e.g., Pezzulli et al., 2005).
% The trend component is a filtered version of the monthly time series.
% 5. The slope of the trend component was calculated using a robust method
% (Sen 1968). The method also calculates the 95% confidence range in the slope.

% The monthly mean values (1) are already tabled.
% Steps 2 & 3, the climatologies and monthly anomalies, is similar to the
% anomalies calculated above. Repeat this, using the 1993-2014 period, to 
% harmonise my outputs with the existing trend map.
% Step 4 is the main difference -- time series method to account for
% seasonality.
% Step 5 is similar to the linear modelling above -- estimating the
% regression coefficient.

% Work with the coarsely resolved data
dat = map_pH_table_9x3(:,1:7);
% Index grid cell by a single column
d = unique(dat(:,1:4));
ngridcells = height(d);
d.gridCell = (1:ngridcells)';
dat = join(dat,d);
% Sort data
[~, I] = sortrows([dat.gridCell, dat.year, dat.month]);
dat = dat(I,:);
% Calculate the climatologies as monthly means over 1993-2014
climatology = nan(12, ngridcells);
climWindow = [1993, 2014];
ind = climWindow(1) <= dat.year & dat.year <= climWindow(2);
for i = 1:ngridcells
    d = dat(ind & dat.gridCell == i,:);
    for j = 1:12
        dm = d(d.month == j,:);
        climatology(j,i) = mean(dm.pH, 'omitnan');
    end
end
% Calculate anomalies as difference from climatology
climatology3d = repmat(reshape(climatology, 12, 1, ngridcells), [1, nyrs, 1]);
dat.anomaly = dat.pH - climatology3d(:);
% Time series analysis using X11 algorithm
dat.seasonal_factor = nan(height(dat), 1);
dat.random_factor = nan(height(dat), 1);
dat.trend = nan(height(dat), 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    vx = x11(d.anomaly, 12, 'add');
    dat.seasonal_factor(ind) = vx.d10;
    dat.random_factor(ind) = vx.d13;
    dat.trend(ind) = vx.d12;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['decompose 9x3 pH data ' num2str(percent_prog) '%'])
    end
end

% Plot time series decomposition outputs
switch plot_time_series_decomposition, case true
    gc = 450; % choose grid cell
    ind = dat.gridCell == gc;
    subplot(2,2,1)
    plot(dat.anomaly(ind))
    title('pH anomaly')
    subplot(2,2,2)
    plot(dat.seasonal_factor(ind))
    title('seasonal factor')
    subplot(2,2,3)
    plot(dat.random_factor(ind))
    title('random factor')
    subplot(2,2,4)
    plot(dat.trend(ind))
    title('trend')
    sgtitle(['grid cell ' num2str(gc)])
end
% Calculate the linear slope of the trend.
% This will produce a single value per grid cell, unlike my previous method
% of calculating separate values for each month and grid cell.
rateOfChange = nan(ngridcells, 1);
pvals = nan(ngridcells, 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    y = d.trend; % temperature anomaly corrected for seasonality
    if all(isnan(y)), continue; end
    x = (1:length(y))' / 12; % independent variable: time (years)
    model = fitlm(x, y); % fit a linear model
    coef = model.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pval = model.Coefficients.pValue(2); % get p-values to gauge significance
    rateOfChange(i) = coef;
    pvals(i) = pval;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 9x3 pH data ' num2str(percent_prog) '%'])
    end
end

d = unique(dat(:,{'lonmin' 'lonmax', 'latmin' 'latmax' 'gridCell'}));
[~,I] = sort(d.gridCell);
d = d(I,:);
d.linearTrend = rateOfChange;
d.pValue = pvals;
dat_trend_9x3 = d;

% Now repeat for higer res data
dat = map_pH_table_3x1(:,1:7);
% Index grid cell by a single column
d = unique(dat(:,1:4));
ngridcells = height(d);
d.gridCell = (1:ngridcells)';
dat = join(dat,d);
% Sort data
[~, I] = sortrows([dat.gridCell, dat.year, dat.month]);
dat = dat(I,:);
% Calculate the climatologies as monthly means over 1993-2014
climatology = nan(12, ngridcells);
climWindow = [1993, 2014];
ind = climWindow(1) <= dat.year & dat.year <= climWindow(2);
for i = 1:ngridcells
    d = dat(ind & dat.gridCell == i,:);
    for j = 1:12
        dm = d(d.month == j,:);
        climatology(j,i) = mean(dm.pH, 'omitnan');
    end
    disp([num2str(100 * i / ngridcells) '%'])
end
% Calculate anomalies as difference from climatology
climatology3d = repmat(reshape(climatology, 12, 1, ngridcells), [1, nyrs, 1]);
dat.anomaly = dat.pH - climatology3d(:);
% Time series analysis using X11 algorithm
dat.seasonal_factor = nan(height(dat), 1);
dat.random_factor = nan(height(dat), 1);
dat.trend = nan(height(dat), 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    vx = x11(d.anomaly, 12, 'add');
    dat.seasonal_factor(ind) = vx.d10;
    dat.random_factor(ind) = vx.d13;
    dat.trend(ind) = vx.d12;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['decompose 3x1 pH data ' num2str(percent_prog) '%'])
    end
end
% Plot time series decomposition outputs
switch plot_time_series_decomposition, case true
    gc = 1000; % choose grid cell
    ind = dat.gridCell == gc;
    subplot(2,2,1)
    plot(dat.anomaly(ind))
    title('temperature anomaly')
    subplot(2,2,2)
    plot(dat.seasonal_factor(ind))
    title('seasonal factor')
    subplot(2,2,3)
    plot(dat.random_factor(ind))
    title('random factor')
    subplot(2,2,4)
    plot(dat.trend(ind))
    title('trend')
    sgtitle(['grid cell ' num2str(gc)])
end
% Calculate the linear slope of the trend.
% This will produce a single value per grid cell, unlike my previous method
% of calculating separate values for each month and grid cell.
rateOfChange = nan(ngridcells, 1);
pvals = nan(ngridcells, 1);
for i = 1:ngridcells
    ind = dat.gridCell == i;
    d = dat(ind,:);
    y = d.trend; % temperature anomaly corrected for seasonality
    if all(isnan(y)), continue; end
    x = (1:length(y))' / 12; % independent variable: time (years)
    model = fitlm(x, y); % fit a linear model
    coef = model.Coefficients.Estimate(2); % extract the trend (degree C / year)
    pval = model.Coefficients.pValue(2); % get p-values to gauge significance
    rateOfChange(i) = coef;
    pvals(i) = pval;
    percent_prog = floor(100 * i / ngridcells);
    percent_prog_ = floor(100 * (i-1) / ngridcells);
    if percent_prog ~= percent_prog_
        disp(['lin. trend 3x1 pH data ' num2str(percent_prog) '%'])
    end
end

d = unique(dat(:,{'lonmin' 'lonmax', 'latmin' 'latmax' 'gridCell'}));
[~,I] = sort(d.gridCell);
d = d(I,:);
d.linearTrend = rateOfChange;
d.pValue = pvals;
dat_trend_3x1 = d;

dat_trend_9x3.gridCell = [];
dat_trend_3x1.gridCell = [];

switch saveOutput, case true
    path = fullfile(dirData, 'pH', 'SOCAT', 'compiled');
    if ~exist(path, 'dir')
        mkdir(path)
            addpath(genpath(path))
    end

    filename = 'pH_time_series_trend_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(dat_trend_9x3, filepath)

    filename = 'pH_time_series_trend_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(dat_trend_3x1, filepath)
end



