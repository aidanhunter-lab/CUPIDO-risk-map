function dat = loadWaterDensity(dataPath, season, months)

% Function that loads water density data from physical model outputs. This
% requires that the density is already calculated for the GLORYS model.

% The data are very large (for whole Southern Ocean) so load, at most, one
% season at a time. A season is defined as starting July 1st and ending June 30th.

Months = ["01","02","03","04","05","06","07","08","09","10","11","12"];
Days = ["15","16"];

nseasons = length(season);
ymin = season;
ymax = season + 1;
yrs = ymin:ymax;
nmonths = length(months);
ntimes = nseasons * nmonths;

% Data are stored with file names given by date, so define a struct with
% the required dates
for i = 1:ntimes
    imonth = mod(i-1, nmonths) + 1;
    month = months(imonth);
    if ismember(month, 7:12)
        iyr = 1;
    elseif ismember(month, 1:6)
        iyr = 2;
    end
    yr = yrs(iyr);
    timeDat.season(i) = season;
    timeDat.year(i) = yr;
    timeDat.month(i) = month;
end

% Load water density data into struct 'dat'
for i = 1:ntimes
    f = {fullfile(dataPath, append(num2str(timeDat.year(i)), Months(timeDat.month(i) == str2double(Months)), Days(1), ".nc")), ...
        fullfile(dataPath, append(num2str(timeDat.year(i)), Months(timeDat.month(i) == str2double(Months)), Days(2), ".nc"))};
    fileExists = [exist(f{1}, 'file'), exist(f{2}, 'file')] == 2;
    if ~any(fileExists), continue; end
    filename = f{fileExists}; % choose the file name that exists in the search path
    if i ==1
%         ncdisp(filename)
        fileInfo = ncinfo(filename);
        % Get data dimensions
        dataDim = [{fileInfo.Dimensions.Name};
            {fileInfo.Dimensions.Length}];
        nlon = dataDim{2,strcmp(dataDim(1,:), 'longitude')};
        nlat = dataDim{2,strcmp(dataDim(1,:), 'latitude')};
        ndep = dataDim{2,strcmp(dataDim(1,:), 'depth')};
        dat.time = nan(ntimes, 1, 'single');
        dat.lon = ncread(filename, 'longitude');
        dat.lat = ncread(filename, 'latitude');
        dat.depth = ncread(filename, 'depth');
        dat.density = nan(nlon, nlat, ndep, ntimes, 'single');
    end
    dat.time(i) = ncread(filename, 'time');
    dat.density(:,:,:,i) = ncread(filename, 'density');
end

% Convert the time variable to MatLab standard time
GLORYS_startDate = [1950 1 1 0 0 0];
adjustTime = datenum(GLORYS_startDate);
dat.time = dat.time / 24; % days since 1/1/1950
dat.time = dat.time + adjustTime;
