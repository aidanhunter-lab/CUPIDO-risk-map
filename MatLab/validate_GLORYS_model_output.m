% Validate GLORYS physical model data using CTD samples

%% Preamble

% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('validate_GLORYS_model_output.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))


%% Load & clean CTD data

% Set data source and region -- in accordance with directory structure
source = 'BODC';
region = 'South Georgia';

[Data, MetaData, ~] = Prepare_CTD_Data(project, source, region, ... 
    'plotMaps', true, 'displayData', false, ...
    'flagGoodMeasure', [49 50 56], 'plotSpuriousPoints', false, ... 
    'removeSpuriousPoints', true, 'onlyDepthProfiles', true, ...
    'ShowMapTitle', true);

% Rename the Label fields -- this could be put into the above function...
cruises = unique(Data.Cruise, 'stable');

for i = 1:length(cruises)
    cruise = cruises{i};
    icruise = contains(MetaData.Cruise, cruise);
    labels = MetaData.Label(icruise);
    ncasts = sum(icruise);
    for j = 1:ncasts
        labels{j} = ['cast' num2str(j)];
    end
    MetaData.Label(icruise) = labels;
    idcruise = contains(Data.Cruise, cruise);
    for j = 1:ncasts
        stationID = MetaData.StationID{icruise & strcmp(MetaData.Label, labels{j})};
        k = idcruise & strcmp(Data.StationID, stationID);
        Data.Label(k) = labels(j);
    end
    clearvars cruise icruise idcruise labels ncasts stationID i j k
end

% for i = 1:length(cruises)
%     cruise = cruises{i};
%     icruise = contains(MetaData.Cruise, cruise);
%     labels = MetaData.StationID(icruise);
%     ncasts = sum(icruise);
%     for j = 1:ncasts
%         labels{j} = ['station' labels{j}];
%     end
%     MetaData.Label(icruise) = labels;
%     idcruise = contains(Data.Cruise, cruise);
%     for j = 1:ncasts
%         stationID = MetaData.StationID{icruise & strcmp(MetaData.Label, labels{j})};
%         k = idcruise & strcmp(Data.StationID, stationID);
%         Data.Label(k) = labels(j);
%     end
%     clearvars cruise icruise idcruise labels ncasts stationID i j k
% end


%% Identify time and location of CTD samples to compare with GLORYS

% To validate the GLORYS model output let's use CTD samples from cruises
% JR20141115 (2014-15), JR15002 (2015-16), and JR16003 (2016-17).
% Omit the other cruises because they either have few samples or else most
% stations are close to each other and near the coast (which causes model
% issues due to boundary conditions).
cruises = {'JR20141115', 'JR15002', 'JR16003'};

% For each cruise we need to find the range of times and locations sampled,
% then match these to the GLORYS model output. To avoid downloading/storing
% excessivley large files, we only retain model output matching the CTD
% samples.

% The BODC data has time variable corresponding to days since 1/1/-4713,
% which differs from MatLab's default of 1/1/0000.
BODC_startDate = [-4712 1 1 0 0 0];
adjustTime = datenum(BODC_startDate);
% Use MatLab's default serial time
MetaData.Time = MetaData.Time + adjustTime;
Data.Time = Data.Time + adjustTime;

% For each cruise, note the GLORYS output data to extract.
% If the GLORYS data is not already downloaded then extractData may be used
% to inform what data is required.
for i = 1:length(cruises)
    % Subset data by cruise
    d = MetaData(strcmp(MetaData.Cruise, cruises(i)),:);
    timeRange = [min(d.Time) max(d.Time)];
    dateRange = datestr(timeRange);
    lonRange = [min(d.Longitude) max(d.Longitude)];
    latRange = [min(d.Latitude) max(d.Latitude)];
    depthRange = [0 max(d.SeaFloorDepth)];
    extractData.(cruises{i}).dates = dateRange;
    extractData.(cruises{i}).lons = lonRange;
    extractData.(cruises{i}).lats = latRange;
    extractData.(cruises{i}).depths = depthRange;
end

clearvars BODC_startDate adjustTime d timeRange dateRange2 lonRange latRange depthRange i


%% Load GLORYS data

% I have the full GLORYS data set for South Georgia stored on my laptop,
% but as it's a large file it may not be suitable for GitHub (this should
% be fine with LFS... but there's some sort of problem).

filename = 'cmems_mod_glo_phy_my_0.083_P1M-m_1656419648164.nc';

ncinfo(filename)
ncdisp(filename)

dat.time = ncread(filename, 'time'); % hours since 1950/1/1
dat.depth = ncread(filename, 'depth'); % metres
dat.longitude = ncread(filename, 'longitude'); % degrees East
dat.latitude = ncread(filename, 'latitude'); % degrees North
dat.temperature = ncread(filename, 'thetao'); % degrees C, seawater potential temperature, add_offset = 21, scale_factor = 0.00073244
dat.salinity = ncread(filename, 'so'); % 1e-3 (practical salinity unit), add_offset = -0.0015259, scale_factor = 0.0015259

% Adjust the time units to match those used for the CTD data, which is
% MatLab default serial time (days since 1/1/0000).
GLORYS_startDate = [1950 1 1 0 0 0];
adjustTime = datenum(GLORYS_startDate);
dat.time = dat.time / 24; % days since 1/1/1950
dat.time = dat.time + adjustTime; % GLORYS time units should now match CTD data

clearvars GLORYS_startDate adjustTime filename


%% Match GLORYS data to CTD samples

% Index space-time points of GLORYS output matching cruises
for i = 1:length(cruises)
    % times
    m = abs(dat.time - datenum(extractData.(cruises{i}).dates)');
    [~, it] = min(m);
    it = it(1):it(2);
    % lat-lons
    m = abs(dat.longitude - extractData.(cruises{i}).lons);
    [~, ilon] = min(m);
    ilon = ilon(1):ilon(2);
    m = abs(dat.latitude - extractData.(cruises{i}).lats);
    [~, ilat] = min(m);
    ilat = ilat(1):ilat(2);
    % depths
    idep = dat.depth < extractData.(cruises{i}).depths(2);
    idep = find(idep + [0; -diff(idep)])';

    indicesGLORYS.(cruises{i}).times = it;
    indicesGLORYS.(cruises{i}).lons = ilon;
    indicesGLORYS.(cruises{i}).lats = ilat;
    indicesGLORYS.(cruises{i}).depths = idep;
    clearvars i m it ilon ilat idep
end

% Filter the GLORYS output data by cruise
for i = 1:length(cruises)
    indices = indicesGLORYS.(cruises{i});
    dat.(cruises{i}).time = dat.time(indices.times);
    dat.(cruises{i}).depth = dat.depth(indices.depths);
    dat.(cruises{i}).longitude = dat.longitude(indices.lons);
    dat.(cruises{i}).latitude = dat.latitude(indices.lats);
    dat.(cruises{i}).temperature = dat.temperature( ...
        indices.lons, indices.lats, indices.depths, indices.times);
    dat.(cruises{i}).salinity = dat.salinity( ...
        indices.lons, indices.lats, indices.depths, indices.times);
    clear vars indices i
end

% Match each individual CTD cast to GLORYS model output
for i = 1:length(cruises)
    md = MetaData(strcmp(MetaData.Cruise, cruises{i}),:);
    ncasts = size(md, 1);
    for j = 1:ncasts
        % time
        m = abs(dat.(cruises{i}).time - md.Time(j));
        [~,it] = min(m);
        % lat-lon
        m = abs(dat.(cruises{i}).longitude - md.Longitude(j));
        [~,ilon] = min(m);
        m = abs(dat.(cruises{i}).latitude - md.Latitude(j));
        [~,ilat] = min(m);
        % depths
        idep = dat.(cruises{i}).depth < md.SeaFloorDepth(j);
        idep = find(idep + [0; -diff(idep)]);

        indicesGLORYS.(cruises{i}).(['cast' num2str(j)]).time = it;
        indicesGLORYS.(cruises{i}).(['cast' num2str(j)]).lon = ilon;
        indicesGLORYS.(cruises{i}).(['cast' num2str(j)]).lat = ilat;
        indicesGLORYS.(cruises{i}).(['cast' num2str(j)]).depths = idep;
    end
    clearvars i j md ncasts it ilon ilat idep
end

for i = 1:length(cruises)
    cruise = cruises{i};
    ncasts = sum(strcmp(MetaData.Cruise, cruises{i}));
    for j = 1:ncasts
        cast = ['cast' num2str(j)];
        indices = indicesGLORYS.(cruise).(cast);
        dat.(cruise).(cast).time = dat.(cruise).time(indices.time);
        dat.(cruise).(cast).depth = dat.(cruise).depth(indices.depths);
        dat.(cruise).(cast).longitude = dat.(cruise).longitude(indices.lon);
        dat.(cruise).(cast).latitude = dat.(cruise).latitude(indices.lat);
        dat.(cruise).(cast).temperature = squeeze(dat.(cruise).temperature( ...
            indices.lon, indices.lat, indices.depths, indices.time));
        dat.(cruise).(cast).salinity = squeeze(dat.(cruise).salinity( ...
            indices.lon, indices.lat, indices.depths, indices.time));
        removeNaNs = isnan(dat.(cruise).(cast).temperature) | ...
            isnan(dat.(cruise).(cast).salinity);
        dat.(cruise).(cast).depth = dat.(cruise).(cast).depth(~removeNaNs);
        dat.(cruise).(cast).temperature = dat.(cruise).(cast).temperature(~removeNaNs);
        dat.(cruise).(cast).salinity = dat.(cruise).(cast).salinity(~removeNaNs);
    end
    clearvars i j cruise ncasts cast indices removeNaNs
end


%% Compare CTD samples to GLORYS model -- make plots

for i = 1:length(cruises)
    cruise = cruises{i};
    ctd = Data(strcmp(Data.Cruise, cruise),:);
    model = dat.(cruise);
    casts = unique(ctd.Label);
    ncasts = length(casts);
    if ncasts <= 6
        ncols = ncasts; nrows = 2;
    else
        ncols = ceil(sqrt(ncasts));
        nrows = 2 * ceil(ncasts / ncols);
    end
    plotHandles.(cruise) = figure;
    set(plotHandles.(cruise), {'Units', 'Position'}, {'inches', [0, 0, 3*ncols, 3*nrows]})
    for j = 1:ncasts
        cast = casts{j};
        station = str2double(MetaData.StationID{...
            strcmp(MetaData.Cruise, cruise) & strcmp(MetaData.Label, cast)});
        ctdj = ctd(strcmp(ctd.Label, cast),:);
        modelj = model.(cast);
        subplot(nrows, ncols, station)
        plot(ctdj.Temperature, ctdj.Depth, 'r');
        hold on
        plot(modelj.temperature, modelj.depth, '--r')
        hold off
        set(gca, 'YDir','reverse')
        xlabel(['temperature (' char(176) 'C)'])
        ylabel('depth (m)')
        title(['station ' num2str(station)], 'FontWeight', 'normal')
    end
    for j = 1:ncasts
        cast = casts{j};
        station = str2double(MetaData.StationID{...
            strcmp(MetaData.Cruise, cruise) & strcmp(MetaData.Label, cast)});
        ctdj = ctd(strcmp(ctd.Label, cast),:);
        modelj = model.(cast);
        subplot(nrows, ncols, station + 0.5 * nrows * ncols)
        plot(ctdj.SalinityPractical, ctdj.Depth, 'b');
        hold on
        plot(modelj.salinity, modelj.depth, '--b')
        hold off
        set(gca, 'YDir','reverse')
        xlabel('salinity')
        ylabel('depth (m)')
        title(['station ' num2str(station)], 'FontWeight', 'normal')
    end
    sgtitle(['Cruise ' cruise])
end


















