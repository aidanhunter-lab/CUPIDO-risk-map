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
sources = {'BODC', 'BAS central storage'};
source = sources{2};

% region = 'South Georgia';
region = [];

switch source
    case 'BODC'
        [Data, MetaData, ~] = Prepare_CTD_Data(project, source, region, ...
            'plotMaps', true, 'displayData', false, ...
            'flagGoodMeasure', [49 50 56], 'plotSpuriousPoints', false, ...
            'removeSpuriousPoints', true, 'onlyDepthProfiles', true, ...
            'ShowMapTitle', true);
    case 'BAS central storage'
        [Data, MetaData, ~] = Prepare_CTD_Data(project, source, region, ...
            'plotMaps', true, 'ShowMapTitle', true);
end


% Rename the Label fields -- this could be put into the above function...
switch source
    case 'BODC'
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
end

%% Identify time and location of CTD samples to compare with GLORYS

switch source
    case 'BODC'
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

        clearvars BODC_startDate adjustTime d timeRange dateRange lonRange latRange depthRange i

    case 'BAS central storage'

        % Find the range of times and locations sampled by CTDs then match
        % these to the GLORYS model output. To avoid downloading/storing
        % excessivley large files, we only retain model output matching the
        % CTD samples.
        timeRange = [min(Data.datenum) max(Data.datenum)];
        dateRange = datestr(timeRange);
        lonRange = [min(Data.lon) max(Data.lon)];
        latRange = [min(Data.lat) max(Data.lat)];
        %         depthRange = Data.press;
        extractData.dates = dateRange;
        extractData.lons = lonRange;
        extractData.lats = latRange;

        clearvars timeRange dateRange lonRange latRange

end



%% Load GLORYS data

% Use info stored in extractData to inform GLORYS data download. Get data
% from: https://resources.marine.copernicus.eu/product-detail/GLOBAL_MULTIYEAR_PHY_001_030/DATA-ACCESS



% I have the full GLORYS data set for South Georgia stored on my laptop,
% but as it's a large file it may not be suitable for GitHub (this should
% be fine with LFS... but there's some sort of problem).

switch source
    case 'BODC'
        filename = 'cmems_mod_glo_phy_my_0.083_P1M-m_1656419648164.nc'; % This South Georgia monthly means, lots of years
    case 'BAS central storage'
        filename = 'cmems_mod_glo_phy_my_0.083_P1D-m_1657890500777.nc'; % Scotia Sea daily means, 2003 Jan-Feb
end

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

switch source
    case 'BODC'
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

    case 'BAS central storage'
        % Index space-time points of GLORYS output to match CTD samples
        ncasts = size(Data.stn, 1);
        for i = 1:ncasts
            % times
            m = abs(dat.time - Data.datenum(i));
            [~,it] = min(m);
            % lat-lon
            m = abs(dat.longitude - Data.lon(i));
            [~,ilon] = min(m);
            m = abs(dat.latitude - Data.lat(i));
            [~,ilat] = min(m);
            % depths
            % I'M NOT SURE HOW TO HANDLE THE DEPTHS BECAUSE THE CTD DATA IS
            % GIVEN IN TERMS OF PRESSURE...
%             idep = dat.depth < Data.press(i);
%             idep = find(idep + [0; -diff(idep)]);
            indicesGLORYS.(['cast' num2str(i)]).time = it;
            indicesGLORYS.(['cast' num2str(i)]).lon = ilon;
            indicesGLORYS.(['cast' num2str(i)]).lat = ilat;
%             indicesGLORYS.(['cast' num2str(i)]).depths = idep;

        end
        clearvars i m it ilon ilat

        for i = 1:ncasts
            cast = ['cast' num2str(i)];
            indices = indicesGLORYS.(cast);
            dat.(cast).time = dat.time(indices.time);
%             dat.(cast).depth = dat.depth(indices.depths);
            dat.(cast).depth = dat.depth;
            dat.(cast).longitude = dat.longitude(indices.lon);
            dat.(cast).latitude = dat.latitude(indices.lat);
            dat.(cast).temperature = squeeze(dat.temperature( ...
                indices.lon, indices.lat, :, indices.time));
            dat.(cast).salinity = squeeze(dat.salinity( ...
                indices.lon, indices.lat, :, indices.time));
            removeNaNs = isnan(dat.(cast).temperature) | isnan(dat.(cast).salinity);
            dat.(cast).depth = dat.(cast).depth(~removeNaNs);
            dat.(cast).temperature = dat.(cast).temperature(~removeNaNs);
            dat.(cast).salinity = dat.(cast).salinity(~removeNaNs);
        end
        clearvars i cast indices removeNaNs

end


%% Compare CTD samples to GLORYS model -- make plots

% Simple plots of measurements and model outputs against depth to visually
% compare CTD measures to physical model.

switch source
    case 'BODC'

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

    case 'BAS central storage'

        nrows = 2; % top = temperature, bottom = salinity
        ncols = 4; % one column per CTD cast
        h = 5; % plot height
        w = 5; % and width
        logDepth = false; % depth axis on log-scale?

        for i = 1:ncasts
            if i == 1, figNum = 0; end
            plotCount = mod(i, ncols);
            if plotCount == 0, plotCount = ncols; end
            newFig = plotCount == 1;
            if newFig
                figNum = figNum + 1;
                figName = ['plot' num2str(figNum)];
                simplePlotHandles.(figName) = figure;
                set(simplePlotHandles.(figName), {'Units', 'Position'}, {'inches', [0, 0, w*ncols, h*nrows]})
            end
            % temperature
            subplot(nrows, ncols, plotCount)
            idat = dat.(['cast' num2str(i)]); % model output corresponding to CTD cast i
            plot(Data.temp(:,i), Data.depth(:,i), 'r');
            hold on
            plot(idat.temperature, idat.depth, '--r')
            hold off
            set(gca, 'YDir','reverse')
            if logDepth
                set(gca, 'YScale', 'log')
            end
            xlabel(['temperature (' char(176) 'C)'])
            ylabel('depth (m)')
            title(['event ' num2str(Data.stn(i))], 'FontWeight', 'normal')
            % salinity
            subplot(nrows, ncols, plotCount + ncols)
            plot(Data.salin(:,i), Data.depth(:,i), 'b');
            hold on
            plot(idat.salinity, idat.depth, '--b')
            hold off
            set(gca, 'YDir','reverse')
            if logDepth
                set(gca, 'YScale', 'log')
            end
            xlabel('salinity')
            ylabel('depth (m)')
            title(['station ' num2str(Data.stn(i))], 'FontWeight', 'normal')
            if plotCount == ncols
                sgtitle('Compare GLORYS model to CTD measurements')
            end
        end
end


close all

% Plot some summary statistics quantifying the GLORYS model 'error' with
% referenced to the CTD measures. For now just plot differences and ratios

% The ratio is an unstable metric because the temperature measurements are
% frequently near zero. This causes large oscillations.
% Use a difference metric, but scaled by the data means. Scaling by the
% data mean will allow the metrics to be comparable across different data
% types, i.e., salinity and temperature.
plotType = 'differences'; % may be 'differences'
% plotType = 'ratios'; % or 'ratios'

comparisons.difTemp = nan(size(Data.temp));
comparisons.difSalin = nan(size(Data.salin));

for i = 1:ncasts
    idat = dat.(['cast' num2str(i)]); % model output corresponding to CTD cast i
    % interpolate model output to match ctd measurement depths
    idat.tempInterp = interp1(idat.depth, idat.temperature, Data.depth(:,i));
    idat.salinInterp = interp1(idat.depth, idat.salinity, Data.depth(:,i));
    comparisons.difTemp(:,i) = idat.tempInterp - Data.temp(:,i);
    comparisons.difSalin(:,i) = idat.salinInterp - Data.salin(:,i);
end


% Find the across-ctd casts mean and st.dev.
summaryStats.mt = mean(comparisons.difTemp, 2, 'omitnan');
summaryStats.ms = mean(comparisons.difSalin, 2, 'omitnan');
summaryStats.sdt = std(comparisons.difTemp, 1, 2, 'omitnan');
summaryStats.sds = std(comparisons.difSalin, 1, 2, 'omitnan');

% Summary plots showing deviation of GLORYS model outputs from CTD measures.
% Temperature

savePlots = true;
lw = 1; % line width for summary stats

fig = figure;
set(fig, {'Units', 'Position'}, {'inches', [0, 0, 8, 4]})
subplot(1,2,1)
Col = [0.65, 0.65, 0.65, 0.3]; % light grey semi-transparent lines for all individual CTD casts
plot(comparisons.difTemp(1:end-1,:), depths, 'Color', Col); % omit the deepest values as they're only 2 measures there so the summary stats are useless
hold on
plot(summaryStats.mt(1:end-1), depths, 'Color', [0, 0, 0], 'LineWidth', lw)
plot([summaryStats.mt(1:end-1) - summaryStats.sdt(1:end-1), summaryStats.mt(1:end-1) + summaryStats.sdt(1:end-1)], ...
    depths, '--', 'Color', [0, 0, 0], 'LineWidth', lw)
xl = xlim;
yl = [0, max(depths)];
plot([0, 0], yl, 'r:', 'LineWidth', lw)
set(gca, 'YDir', 'reverse')
ylabel('depth (m)')
xlabel('temperature difference')
yl = ylim;

l = 1/8*diff(xl); % legend line length
lx = 1/12*diff(xl); % legend x inset
ly = 1/20*diff(yl); % legend y inset
lys = 1/20*diff(yl); % legend y spacing
lts = 0.25 * l; % text space
fs = 8; % font size

plot([xl(1)+lx, xl(1)+lx+l], [yl(2)-ly, yl(2)-ly], 'r:', 'LineWidth', lw)
text(lts + xl(1)+lx+l, yl(2)-ly, 'zero difference', 'FontSize', fs)
plot([xl(1)+lx, xl(1)+lx+l], [yl(2)-ly-lys, yl(2)-ly-lys], 'Color', Col)
text(lts + xl(1)+lx+l, yl(2)-ly-lys, 'individual samples', 'FontSize', fs)
plot([xl(1)+lx, xl(1)+lx+l], [yl(2)-ly-2*lys, yl(2)-ly-2*lys], '--', 'Color', [0, 0, 0], 'LineWidth', lw)
text(lts + xl(1)+lx+l, yl(2)-ly-2*lys, 'mean \pm st.dev.', 'FontSize', fs)
plot([xl(1)+lx, xl(1)+lx+l], [yl(2)-ly-3*lys, yl(2)-ly-3*lys], 'Color', [0, 0, 0], 'LineWidth', lw)
text(lts + xl(1)+lx+l, yl(2)-ly-3*lys, 'mean', 'FontSize', fs)
hold off

% Salinity
subplot(1,2,2)
Col = [0.65, 0.65, 0.65, 0.3]; % light grey semi-transparent lines for all individual CTD casts
plot(comparisons.difSalin(1:end-1,:), depths, 'Color', Col); % omit the deepest values as they're only 2 measures there so the summary stats are useless
hold on
plot(summaryStats.ms(1:end-1), depths, 'Color', [0, 0, 0], 'LineWidth', lw)
plot([summaryStats.ms(1:end-1) - summaryStats.sds(1:end-1), summaryStats.ms(1:end-1) + summaryStats.sds(1:end-1)], ...
    depths, '--', 'Color', [0, 0, 0], 'LineWidth', lw)
yl = [0, max(depths)];
% yl = ylim;
plot([0, 0], yl, 'r:', 'LineWidth', lw)
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'XLim', [-0.5, 0.5]) % crop the plot
xlabel('salinity difference')

sgtitle({'Compare GLORYS model output to CTD measures', 'Lines are modelled minus measured values'})

switch savePlots, case true
    filename = 'GLORYS vs CTD JR82_salinity and temperature.png';
    filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
    exportgraphics(fig, filepath)
end





switch source
    case 'BAS central storage'

        nrows = 2; % top = temperature, bottom = salinity
        ncols = 4; % one column per CTD cast
        h = 5; % plot height
        w = 5; % and width
        logDepth = false; % depth axis on log-scale?

        for i = 1:ncasts
            if i == 1, figNum = 0; end
            plotCount = mod(i, ncols);
            if plotCount == 0, plotCount = ncols; end
            newFig = plotCount == 1;
            if newFig
                figNum = figNum + 1;
                figName = ['plot' num2str(figNum)];
                diffPlotHandles.(figName) = figure;
                set(diffPlotHandles.(figName), {'Units', 'Position'}, {'inches', [0, 0, w*ncols, h*nrows]})
            end
            % temperature
            subplot(nrows, ncols, plotCount)
            idat = dat.(['cast' num2str(i)]); % model output corresponding to CTD cast i
            dataTemp = Data.temp(:,i);
            modTemp = idat.temperature;
            dataDepth = Data.depth(:,i);
            modDepth = idat.depth;
            % interpolate model output to match ctd measurement depths
            modTempi = interp1(modDepth, modTemp, dataDepth);
            difTemp = modTempi - dataTemp;
            ratTemp = modTempi ./ dataTemp;
            switch plotType
                case 'ratios'
                    xvar = ratTemp;
                case 'differences'
                    xvar = difTemp;
%                     xvar = xvar ./ mean(dataTemp, 'omitnan');
%                     xvar = xvar ./ std(dataTemp, 'omitnan');
            end
            plot(xvar, dataDepth, 'r')
            set(gca, 'YDir','reverse')
            if logDepth
                set(gca, 'YScale', 'log')
            end
            xlabel(['temperature difference (' char(176) 'C)'])
            ylabel('depth (m)')
            title(['event ' num2str(Data.stn(i))], 'FontWeight', 'normal')
            % salinity
            subplot(nrows, ncols, plotCount + ncols)
            dataSal = Data.salin(:,i);
            modSal = idat.salinity;
            dataDepth = Data.depth(:,i);
            modDepth = idat.depth;
            % interpolate model output to match ctd measurement depths
            modSali = interp1(modDepth, modSal, dataDepth);
            difSal = modSali - dataSal;
            ratSal = modSali ./ dataSal;
            switch plotType
                case 'ratios'
                    xvar = ratSal;
                case 'differences'
                    xvar = difSal;
%                     xvar = xvar ./ mean(dataSal, 'omitnan');
%                     xvar = xvar ./ std(dataSal, 'omitnan');
            end
            plot(xvar, dataDepth, 'b')
            set(gca, 'YDir','reverse')
            if logDepth
                set(gca, 'YScale', 'log')
            end
            xlabel('salinity difference')
            ylabel('depth (m)')
            title(['station ' num2str(Data.stn(i))], 'FontWeight', 'normal')
            if plotCount == ncols
                sgtitle('Compare GLORYS model to CTD measurements')
            end
        end
end

