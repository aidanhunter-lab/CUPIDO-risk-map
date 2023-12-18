%% Download satellite SST data directly from MatLab

% This requires Python is installed
pythonVersion = '3';

% For instructions and dependencies see:
% https://help.marine.copernicus.eu/en/articles/4799385-can-i-download-copernicus-marine-data-via-r-or-matlab#h_db0930094c

%% Copernicus Marine Credentials

% Store your Copernicus username & password in a hidden text file called
% '.Copernicus Marine Credentials.txt', formatted as a table with column
% names 'user' & 'pwd'.
cmc = readtable('.Copernicus Marine Credentials.txt');

% Insert your username
user = cmc.user {1};
% and your password
pwd = cmc.pwd{1};

%% Download parameters for ESA SST CCI and C3S reprocessed sea surface temperature analyses

% Download server
motu_server = 'https://my.cmems-du.eu/motu-web/Motu';
% Model services
serviceID = 'SST_GLO_SST_L4_REP_OBSERVATIONS_010_024-TDS';
% Model product
productID = 'ESACCI-GLO-SST-L4-REP-OBS-SST';


%% Directories
thisFile = which('download_SST.m');
dirBase = fileparts(fileparts(thisFile));
dirData = fullfile(dirBase, 'data', 'sst', 'ESA');
addpath(genpath(dirData))

out_dir = fullfile(dirData);
% If directory does not exist then create it
dir_exists = exist(out_dir, 'dir') == 7;
if ~dir_exists, mkdir(dirData); end

% Code spaces with a backslash (Linux style)
out_dir = strrep(out_dir, ' ', '\ ');


%% Choose download variables
variables = "analysed_sst";

%% Download the data

% Build the command line for motuclient.
% Output variables
vars_ = arrayfun(@(z) sprintf('%s ', '--variable', z), variables, 'UniformOutput', false);
vars = ' ';
for i = 1:length(variables), vars = append(vars, vars_{i}); end
vars = vars(1:end-1);
% Coordinates
lonrange = [-180, 180];
latrange = [-90, -45];
lonrange_ = string(lonrange);
latrange_ = string(latrange);

% Dates - these are daily data spanning a few decades, with a high spatial
% resolution. Download the daily files then takes averages before saving
% them. I need to choose a spatiotemporal resolution. Monthly means and 3*1
% degree lon*lat grid should be OK.
reslon_coarse = 9; % coarse resolution matches krill map
reslat_coarse = 3;
reslon_high = 3; % high (3*1) res is good for interactive map
reslat_high = 1;
longrid_coarse = lonrange(1):reslon_coarse:lonrange(2);
latgrid_coarse = latrange(1):reslat_coarse:latrange(2);
longrid_high = lonrange(1):reslon_high:lonrange(2);
latgrid_high = latrange(1):reslat_high:latrange(2);
nlon_coarse = -1 + length(longrid_coarse);
nlat_coarse = -1 + length(latgrid_coarse);
nlon_high = -1 + length(longrid_high);
nlat_high = -1 + length(latgrid_high);

yrrange = [1982, 2016];
startdate = datetime([yrrange(1) 1 1]);
enddate = datetime([yrrange(2) 12 31]);
alldays = datenum(startdate):datenum(enddate);
% ndays = length(alldays);
[yr, mon, day] = datevec(alldays);
% group data by month
dm = [0, diff(mon)]; dm(dm == -11) = 1;
group = 1 + cumsum(dm);
yr = yr(:); mon = mon(:); day = day(:); group = group(:);
alldates = table(yr, mon, day, group);
alldates.day1 = cell(height(alldates),1); alldates.day1(:) = {'01'};
alldates.day2 = cell(height(alldates),1);
ngroups = max(alldates.group);
for i = 1:ngroups
    j = alldates.group == i;
    alldates.day2(j) = {num2str(max(alldates.day(j)))};
end

alldates.day = [];
alldates = unique(alldates, 'stable');


% Generate & store the Python code to download the daily data
ml1 = append("python", pythonVersion, " -m motuclient --motu ", ...
        motu_server, " --service-id ", serviceID, " --product-id ", productID, ...
        " --longitude-min ", lonrange_(1), " --longitude-max ", lonrange_(2), ...
        " --latitude-min ", latrange_(1), " --latitude-max ", latrange_(2));
ml2 = append(vars, " --out-dir ", out_dir, " --out-name ");
ml3 = append(" --user ", user, " --pwd ", pwd);
mon0 = alldates.mon < 10;
mon = string(alldates.mon);
mon(mon0) = arrayfun(@(z) append('0', z), mon(mon0));
% day0 = alldates.day < 10;
% day = string(alldates.day);
% day(day0) = arrayfun(@(z) append('0', z), day(day0));
yr = string(alldates.yr);
date = [append('"', yr, '-', mon, '-', alldates.day1, " 00:00:00", '"'), ...
    append('"', yr, '-', mon, '-', alldates.day2, " 23:59:59", '"')];
alldates.filename = append(yr, mon, '.nc');
alldates.motu = append(ml1, " --date-min ", date(:,1), " --date-max ", date(:,2), ...
    ml2, alldates.filename, ml3);

% Download data for each month separately, take averages, save the data we
% need and discard the rest, then repeat for each month.
blank_coarse = nan(nlon_coarse, nlat_coarse);
blank_high = nan(nlon_high, nlat_high);
SST_coarse = nan(nlon_coarse, nlat_coarse, ngroups);
SST_high = nan(nlon_high, nlat_high, ngroups);
tic
warn = false(1, ngroups);
for i = 1:ngroups
    if i == 1, wb = waitbar(0, 'Progress...'); end
    d = alldates(i,:);
    % Download data for one month
    code = d.motu;
    s = system(code);
    % store iterations where download failed
    if ~exist('s', 'var'), warn(i) = true; continue; end
    if s ~= 0, warn(i) = true; continue; end
    clearvars s
    % Load data into workspace
    filename = d.filename;
    filepath = fullfile(out_dir, filename);
%     ncinfo(filepath)
%     ncdisp(filepath)
    lon = ncread(filepath, 'lon');
    lat = ncread(filepath, 'lat');
%     time = ncread(filepath, 'time');
    sst = ncread(filepath, 'analysed_sst');
    % Average the daily values
    sst = mean(sst, 3);
    % Average the high-res data into coarser grids
    sst_coarse = blank_coarse;
    sst_high = blank_high;
    for j = 1:nlon_coarse
        jj = longrid_coarse(j) < lon & lon < longrid_coarse(j+1);
        for k = 1:nlat_coarse
            kk = latgrid_coarse(k) < lat & lat < latgrid_coarse(k+1);
            sst_coarse(j,k) = mean(sst(jj & kk'), 'omitnan');
        end
    end
    for j = 1:nlon_high
        jj = longrid_high(j) < lon & lon < longrid_high(j+1);
        for k = 1:nlat_high
            kk = latgrid_high(k) < lat & lat < latgrid_high(k+1);
            sst_high(j,k) = mean(sst(jj & kk'), 'omitnan');
        end
    end
    SST_coarse(:,:,i) = sst_coarse;
    SST_high(:,:,i) = sst_high;
    clearvars sst_coarse sst_high
    % Delete the high-res data
    delete(filepath)
    waitbar(i/ngroups, wb)
end
toc

% Save the averaged data
lon_coarse = 0.5 .* (longrid_coarse(1:end-1) + longrid_coarse(2:end)); lon_coarse = lon_coarse(:);
lat_coarse = 0.5 .* (latgrid_coarse(1:end-1) + latgrid_coarse(2:end)); lat_coarse = lat_coarse(:);
lon_high = 0.5 .* (longrid_high(1:end-1) + longrid_high(2:end)); lon_high = lon_high(:);
lat_high = 0.5 .* (latgrid_high(1:end-1) + latgrid_high(2:end)); lat_high = lat_high(:);

filepath_coarse = fullfile(out_dir,...
    ['sst_monthly_means_1982_2016_res_' num2str(reslon_coarse) 'x' num2str(reslat_coarse) '.nc']);
filepath_high = fullfile(out_dir,...
    ['sst_monthly_means_1982_2016_res_' num2str(reslon_high) 'x' num2str(reslat_high) '.nc']);

nccreate(filepath_coarse, 'analysed_sst', 'Dimensions', {'lon' nlon_coarse 'lat' nlat_coarse 'time' ngroups})
nccreate(filepath_coarse, 'time', 'Dimensions', {'time' ngroups})
nccreate(filepath_coarse, 'lat', 'Dimensions', {'lat' nlat_coarse})
nccreate(filepath_coarse, 'lon', 'Dimensions', {'lon' nlon_coarse})
ncwrite(filepath_coarse, 'analysed_sst', SST_coarse)
ncwrite(filepath_coarse, 'time', 0:(ngroups-1))
ncwrite(filepath_coarse, 'lat', lat_coarse)
ncwrite(filepath_coarse, 'lon', lon_coarse)

nccreate(filepath_high, 'analysed_sst', 'Dimensions', {'lon' nlon_high 'lat' nlat_high 'time' ngroups})
nccreate(filepath_high, 'time', 'Dimensions', {'time' ngroups})
nccreate(filepath_high, 'lat', 'Dimensions', {'lat' nlat_high})
nccreate(filepath_high, 'lon', 'Dimensions', {'lon' nlon_high})
ncwrite(filepath_high, 'analysed_sst', SST_high)
ncwrite(filepath_high, 'time', 0:(ngroups-1))
ncwrite(filepath_high, 'lat', lat_high)
ncwrite(filepath_high, 'lon', lon_high)

% Clear the large high-res data files from the Git Repo now that they're
% averaged
for i = 1:ngroups
    f = fullfile(dirData, alldates.filename(i));
    if isfile(f)
        delete(f)
    end
end



%% Download more recent data

% python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id SST_GLO_SST_L4_REP_OBSERVATIONS_010_024-TDS --product-id C3S-GLO-SST-L4-REP-OBS-SST --longitude-min -180 --longitude-max 180 --latitude-min -90 --latitude-max -45 --date-min "2018-01-01 00:00:00" --date-max "2018-01-01 23:59:59" --variable analysed_sst --variable analysed_sst_uncertainty --variable mask --variable sea_ice_fraction --out-dir <OUTPUT_DIRECTORY> --out-name <OUTPUT_FILENAME> --user <USERNAME> --pwd <PASSWORD>

% Model product
productID = 'C3S-GLO-SST-L4-REP-OBS-SST';


yrrange = [2017, 2021];
startdate = datetime([yrrange(1) 1 1]);
enddate = datetime([yrrange(2) 12 31]);
alldays = datenum(startdate):datenum(enddate);
ndays = length(alldays);
[yr, mon, day] = datevec(alldays);
% group data by month
dm = [0, diff(mon)]; dm(dm == -11) = 1;
group = 1 + cumsum(dm);
yr = yr(:); mon = mon(:); day = day(:); group = group(:);
alldates = table(yr, mon, day, group);
alldates.day1 = cell(height(alldates),1); alldates.day1(:) = {'01'};
alldates.day2 = cell(height(alldates),1);
ngroups = max(alldates.group);
for i = 1:ngroups
    j = alldates.group == i;
    alldates.day2(j) = {num2str(max(alldates.day(j)))};
end
alldates.day = [];
alldates = unique(alldates, 'stable');

% Generate & store the Python code to download the daily data
ml1 = append("python", pythonVersion, " -m motuclient --motu ", ...
        motu_server, " --service-id ", serviceID, " --product-id ", productID, ...
        " --longitude-min ", lonrange_(1), " --longitude-max ", lonrange_(2), ...
        " --latitude-min ", latrange_(1), " --latitude-max ", latrange_(2));
ml2 = append(vars, " --out-dir ", out_dir, " --out-name ");
ml3 = append(" --user ", user, " --pwd ", pwd);
mon0 = alldates.mon < 10;
mon = string(alldates.mon);
mon(mon0) = arrayfun(@(z) append('0', z), mon(mon0));
yr = string(alldates.yr);
date = [append('"', yr, '-', mon, '-', alldates.day1, " 00:00:00", '"'), ...
    append('"', yr, '-', mon, '-', alldates.day2, " 23:59:59", '"')];
alldates.filename = append(yr, mon, '.nc');
alldates.motu = append(ml1, " --date-min ", date(:,1), " --date-max ", date(:,2), ...
    ml2, alldates.filename, ml3);

% Download data for each month separately, take averages, save the data we
% need and discard the rest, then repeat for each month.
blank_coarse = nan(nlon_coarse, nlat_coarse);
blank_high = nan(nlon_high, nlat_high);
SST_coarse = nan(nlon_coarse, nlat_coarse, ngroups);
SST_high = nan(nlon_high, nlat_high, ngroups);
tic
warn = false(1, ngroups);
for i = 1:ngroups
    if i == 1, wb = waitbar(0, 'Progress...'); end
    d = alldates(i,:);
    % Download data for one month
    code = d.motu;
    s = system(code);
    % store iterations where download failed
    if ~exist('s', 'var'), warn(i) = true; continue; end
    if s ~= 0, warn(i) = true; continue; end
    % Load data into workspace
    filename = d.filename;
    filepath = fullfile(out_dir, filename);
%     ncinfo(filepath)
%     ncdisp(filepath)
    lon = ncread(filepath, 'lon');
    lat = ncread(filepath, 'lat');
%     time = ncread(filepath, 'time');
    sst = ncread(filepath, 'analysed_sst');
    % Average the daily values
    sst = mean(sst, 3);
    % Average the high-res data into coarser grids
    sst_coarse = blank_coarse;
    sst_high = blank_high;
    for j = 1:nlon_coarse
        jj = longrid_coarse(j) < lon & lon < longrid_coarse(j+1);
        for k = 1:nlat_coarse
            kk = latgrid_coarse(k) < lat & lat < latgrid_coarse(k+1);
            sst_coarse(j,k) = mean(sst(jj & kk'), 'omitnan');
        end
    end
    for j = 1:nlon_high
        jj = longrid_high(j) < lon & lon < longrid_high(j+1);
        for k = 1:nlat_high
            kk = latgrid_high(k) < lat & lat < latgrid_high(k+1);
            sst_high(j,k) = mean(sst(jj & kk'), 'omitnan');
        end
    end
    SST_coarse(:,:,i) = sst_coarse;
    SST_high(:,:,i) = sst_high;
    clearvars sst_coarse sst_high
    % Delete the high-res data
    delete(filepath)
    waitbar(i/ngroups, wb)
end
toc

% Save the averaged data
lon_coarse = 0.5 .* (longrid_coarse(1:end-1) + longrid_coarse(2:end)); lon_coarse = lon_coarse(:);
lat_coarse = 0.5 .* (latgrid_coarse(1:end-1) + latgrid_coarse(2:end)); lat_coarse = lat_coarse(:);
lon_high = 0.5 .* (longrid_high(1:end-1) + longrid_high(2:end)); lon_high = lon_high(:);
lat_high = 0.5 .* (latgrid_high(1:end-1) + latgrid_high(2:end)); lat_high = lat_high(:);

filepath_coarse = fullfile(out_dir,...
    ['sst_monthly_means_2017_2021_res_' num2str(reslon_coarse) 'x' num2str(reslat_coarse) '.nc']);
filepath_high = fullfile(out_dir,...
    ['sst_monthly_means_2017_2021_res_' num2str(reslon_high) 'x' num2str(reslat_high) '.nc']);

nccreate(filepath_coarse, 'analysed_sst', 'Dimensions', {'lon' nlon_coarse 'lat' nlat_coarse 'time' ngroups})
nccreate(filepath_coarse, 'time', 'Dimensions', {'time' ngroups})
nccreate(filepath_coarse, 'lat', 'Dimensions', {'lat' nlat_coarse})
nccreate(filepath_coarse, 'lon', 'Dimensions', {'lon' nlon_coarse})
ncwrite(filepath_coarse, 'analysed_sst', SST_coarse)
ncwrite(filepath_coarse, 'time', 0:(ngroups-1))
ncwrite(filepath_coarse, 'lat', lat_coarse)
ncwrite(filepath_coarse, 'lon', lon_coarse)

nccreate(filepath_high, 'analysed_sst', 'Dimensions', {'lon' nlon_high 'lat' nlat_high 'time' ngroups})
nccreate(filepath_high, 'time', 'Dimensions', {'time' ngroups})
nccreate(filepath_high, 'lat', 'Dimensions', {'lat' nlat_high})
nccreate(filepath_high, 'lon', 'Dimensions', {'lon' nlon_high})
ncwrite(filepath_high, 'analysed_sst', SST_high)
ncwrite(filepath_high, 'time', 0:(ngroups-1))
ncwrite(filepath_high, 'lat', lat_high)
ncwrite(filepath_high, 'lon', lon_high)


% Clear the large high-res data files from the Git Repo now that they're
% averaged
for i = 1:ngroups
    f = fullfile(dirData, alldates.filename(i));
    if isfile(f)
        delete(f)
    end
end



%% Combine data from the two time spans: 1982-2016 & 2017-2021

filename1_3x1 = 'sst_monthly_means_1982_2016_res_3x1.nc';
filename2_3x1 = 'sst_monthly_means_2017_2021_res_3x1.nc';
filename1_9x3 = 'sst_monthly_means_1982_2016_res_9x3.nc';
filename2_9x3 = 'sst_monthly_means_2017_2021_res_9x3.nc';
filepath1_3x1 = fullfile(dirData, filename1_3x1);
filepath2_3x1 = fullfile(dirData, filename2_3x1);
filepath1_9x3 = fullfile(dirData, filename1_9x3);
filepath2_9x3 = fullfile(dirData, filename2_9x3);

% Handle the coarsely resolved data first
% lon1 = ncread(filepath1_9x3, 'lon');
% lat1 = ncread(filepath1_9x3, 'lat');
time1 = ncread(filepath1_9x3, 'time');
sst1 = ncread(filepath1_9x3, 'analysed_sst');
% lon2 = ncread(filepath2_9x3, 'lon');
% lat2 = ncread(filepath2_9x3, 'lat');
time2 = ncread(filepath2_9x3, 'time');
sst2 = ncread(filepath2_9x3, 'analysed_sst');

time = [time1; time2 + time1(end) + 1]; % months since Jan 1982
nmonths = length(time);
sst = cat(3, sst1, sst2);

filepath = fullfile(out_dir,...
    ['sst_monthly_means_1982_2021_res_' num2str(reslon_coarse) 'x' num2str(reslat_coarse) '.nc']);

nccreate(filepath, 'analysed_sst', 'Dimensions', {'lon' nlon_coarse 'lat' nlat_coarse 'time' nmonths})
nccreate(filepath, 'time', 'Dimensions', {'time' nmonths})
nccreate(filepath, 'lat', 'Dimensions', {'lat' nlat_coarse})
nccreate(filepath, 'lon', 'Dimensions', {'lon' nlon_coarse})
ncwrite(filepath, 'analysed_sst', sst)
ncwrite(filepath, 'time', time)
ncwrite(filepath, 'lat', lat_coarse)
ncwrite(filepath, 'lon', lon_coarse)

% Now the higher resolution data
% lon1 = ncread(filepath1_3x1, 'lon');
% lat1 = ncread(filepath1_3x1, 'lat');
time1 = ncread(filepath1_3x1, 'time');
sst1 = ncread(filepath1_3x1, 'analysed_sst');
% lon2 = ncread(filepath2_3x1, 'lon');
% lat2 = ncread(filepath2_3x1, 'lat');
time2 = ncread(filepath2_3x1, 'time');
sst2 = ncread(filepath2_3x1, 'analysed_sst');

time = [time1; time2 + time1(end) + 1]; % months since Jan 1982
nmonths = length(time);
sst = cat(3, sst1, sst2);

filepath = fullfile(out_dir,...
    ['sst_monthly_means_1982_2021_res_' num2str(reslon_high) 'x' num2str(reslat_high) '.nc']);

nccreate(filepath, 'analysed_sst', 'Dimensions', {'lon' nlon_high 'lat' nlat_high 'time' nmonths})
nccreate(filepath, 'time', 'Dimensions', {'time' nmonths})
nccreate(filepath, 'lat', 'Dimensions', {'lat' nlat_high})
nccreate(filepath, 'lon', 'Dimensions', {'lon' nlon_high})
ncwrite(filepath, 'analysed_sst', sst)
ncwrite(filepath, 'time', time)
ncwrite(filepath, 'lat', lat_high)
ncwrite(filepath, 'lon', lon_high)


