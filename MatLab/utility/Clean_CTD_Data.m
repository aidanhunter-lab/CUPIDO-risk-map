function [Data, MetaData, cruises] = Clean_CTD_Data(project, source, region, ...
    Data, MetaData, cruises, varargin)
% Filter the data to remove spurious points & organise

extractVarargin(varargin)

if ~exist('plotMaps', 'var'), plotMaps = true; end
if ~exist('flagGoodMeasure', 'var'), flagGoodMeasure = [49 50 56]; end % examine NetCDF files to check measurement quality control flags
if ~exist('plotSpuriousPoints', 'var'), plotSpuriousPoints = true; end
if ~exist('removeSpuriousPoints', 'var'), removeSpuriousPoints = true; end
if ~exist('onlyDepthProfiles', 'var'), onlyDepthProfiles = false; end
% parameter used to identify depth profiles
if ~exist('CoV_lim', 'var'), CoV_lim = 0.05; end
% parameters used to (loosely) identify spurious measurements
if ~exist('cutOffValue', 'var'), cutOffValue = 0.5; end                    % test moving std against cutOffValue
if ~exist('leeway', 'var'), leeway = 2; end                              % non-negative integer -- extra data points to test before specifying a point as spurious.
if ~exist('movstdRange', 'var'), movstdRange = 0.1; end                    % range of moving-std as fraction of total points
if ~exist('XaxisLocation', 'var'), XaxisLocation = 'bottom'; end
if ~exist('ShowMapTitle', 'var'), ShowMapTitle = true; end

switch source
    case 'BODC'
        %% Clean the data: remove NaNs & check quality control flagged measures

        % Measurements include some NaN values -- remove these.
        % Then set to NaN any measurements with quality control warning

        cruises_ = rmfield(cruises, 'ncruises');
        fields = fieldnames(cruises_)';
        seasons.nseasons = length(fields);
        seasons.season = cellfun(@(z) strrep(z(7:end), '_', '-'), fields, 'UniformOutput', false);
        seasons.fieldName = fields;
        clearvars cruises_ fields

        % Omit sample depths associated with NaNs in the key measurements of
        % temperature & salinity -- other variables are also useful, but less
        % important
        checkVars = {'Temperature', 'SalinityPractical', 'Depth'};
        for i = 1:seasons.nseasons
            seasonID = seasons.fieldName{i};
            dc = cruises.(seasonID);
            for j = 1:dc.ncruises
                cruiseID = dc.ID{j};
                Dat = Data.(seasonID).(cruiseID);
                samples = fieldnames(Dat);
                nsamples = length(samples);
                for k = 1:nsamples
                    dat = Dat.(samples{k});
                    nz = length(dat.HeightAboveSeaFloor);
                    fields = fieldnames(dat); nfields = length(fields);
                    keepHeight = true(nz, 1); % index all retained sample depths -- those without NaN measurement values
                    for l = 1:nfields
                        if ~any(contains(checkVars, fields{l})), continue; else
                            x = dat.(fields{l});
                            if ~isnumeric(x), continue; end
                            if size(x, 1) ~= nz, continue; else
                                keepHeight = keepHeight & ~isnan(x);
                            end
                        end
                    end
                    for l = 1:nfields
                        x = dat.(fields{l});
                        if ~isnumeric(x), continue; end
                        if size(x, 1) ~= nz, continue
                        else, dat.(fields{l}) = x(keepHeight); end
                    end
                    Data.(seasonID).(cruiseID).(samples{k}) = dat;
                end
            end
        end
        % Rename CTD samples -- in case any were removed above
        for i = 1:seasons.nseasons
            seasonID = seasons.fieldName{i};
            for j = 1:cruises.(seasonID).ncruises
                cruiseID = cruises.(seasonID).ID{j};
                Dat = Data.(seasonID).(cruiseID);
                samples = fieldnames(Dat);
                for k = 1:length(samples)
                    if ~isempty(Dat.(samples{k}).Depth), continue; end
                    Dat = rmfield(Dat, samples{k});
                end
                samples = fieldnames(Dat);
                Dat_ = Dat; clearvars Dat
                for k = 1:length(samples)
                    CTD = ['CTD' num2str(k)];
                    Dat.(CTD) = Dat_.(samples{k});
                    Dat.(CTD).Label = CTD;
                end
                Data.(seasonID).(cruiseID) = Dat;
                k = strcmp(MetaData.Cruise, cruiseID);
                md = MetaData(k,:);
                md = md(ismember(md.Label, samples),:);
                for l = 1:height(md), md.Label{l} = ['CTD' num2str(l)]; end
                MetaData(k,:) = []; MetaData = [MetaData; md];
            end
        end
        clearvars checkVars Dat Dat_ dat dc md nz fields nfields keepHeight samples nsamples CTD x i j k l *ID


        % Set to NaN any measurements with quality control flags
        % flagGoodMeasure = [49 50 56]; % maybe include 51 because it seems that too often depth is a 'probably bad value'... IS THIS RELATED TO THE APPARENT LACK OF DEPTH PROFILES?
        for i = 1:seasons.nseasons
            seasonID = seasons.fieldName{i};
            dc = cruises.(seasonID);
            for j = 1:dc.ncruises
                cruiseID = dc.ID{j};
                Dat = Data.(seasonID).(cruiseID);
                samples = fieldnames(Dat);
                nsamples = length(samples);
                for k = 1:nsamples
                    sampleID = samples{k};
                    dat = Dat.(sampleID);
                    fields = fieldnames(dat);
                    isflag = contains(fields, 'Flag'); % index quality control flag variables
                    for l = 1:length(fields)
                        if ~isflag(l), continue
                        else
                            f = dat.(fields{l});
                            trusted = ismember(f, flagGoodMeasure);
                            if all(trusted), continue
                            else
                                vn = strrep(fields{l}, 'Flag', ''); % variable name
                                x = dat.(vn);
                                x(~trusted) = nan; % measurements with quality control warning set to NaN
                                dat.(vn) = x;
                            end
                        end
                    end
                    % Omit all measurements where depth is NaN
                    isNaN = isnan(dat.Depth);
                    if any(isNaN)
                        for l = 1:length(fields)
                            if size(dat.(fields{l}), 1) ~= length(isNaN), continue; end
                            dat.(fields{l}) = dat.(fields{l})(~isNaN);
                        end
                    end
                    Data.(seasonID).(cruiseID).(sampleID) = dat;
                end
            end
        end
        clearvars dc Dat dat samples nsamples fields f flagGoodMeasure isflag trusted vn x i j k l *ID

        % Note that spurious measurements remain in most samples... It seems that
        % the quality control flags regard initial CTD measures as OK despite the
        % clearly visible dodgyness of initial readings. See example testPlot below
        % and play around with omitting initial measures... Clearly there's more
        % data clearning to do. I'll come back to this later.

        testPlot = figure;
        d = Data.season2016_2017.JR16003.CTD4; % choose, arbitrarily, a CTD sample to plot
        v = 'Temperature';
        scatter(d.(v), -d.Depth)
        i = 1:51; % index spurious points corresponding to initial readings
        hold on
        scatter(d.Temperature(i), -d.Depth(i), 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0])
        xlabel(v); ylabel('depth (m)')
        ca = gca;
        set(ca, 'YTickLabel', strrep(ca.YTickLabel, '-', ''))
        title({'Example plot: red points are spurious initial CTD readings.', 'These should be removed somehow...'})

        close(testPlot)
        clearvars testPlot d v i ca


        %% View CTD sample locations on map
        coordsTable = readtable('regional bounding coordinates.csv');

        if plotMaps
            % Create a separate map for each cruise
            for i = 1:seasons.nseasons
                if i == 1, counter = 0; end
                seasonID = seasons.fieldName{i};
                for j = 1:cruises.(seasonID).ncruises
                    counter = counter + 1;
                    cruiseID = cruises.(seasonID).ID{j};
                    Dat = Data.(seasonID).(cruiseID);
                    samples = fieldnames(Dat);
                    nsamples = length(samples);
                    % Extract CTD sample locations & create ID tags
                    lons = nan(nsamples, 1); lats = nan(nsamples, 1); labels = cell(nsamples, 1);
                    for k = 1:nsamples
                        sampleID = samples{k};
                        lons(k) = Dat.(sampleID).Longitude;
                        lats(k) = Dat.(sampleID).Latitude;
                        l = Dat.(sampleID).Label;
                        labels{k} = num2str(str2double(l(4:end)));
                    end
                    % Create South Georgia base map
                    mapName = ['map' num2str(counter)];
                    figure
                    plotHandles.(mapName) = gcf;
                    %             assignin('base', mapName, figure)
                    plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                        'XaxisLocation', XaxisLocation);
                    % Extract map bounding coordinates
                    [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0]);
                    mapAx = gca;
                    % Plot sample positions
                    hold(mapAx, 'on')
                    m_scatter(lons, lats, 18, [1 0 0], 'filled')
                    m_text(lons, lats, labels, 'FontSize', 8, ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                    titleLon = mean(lonMap);
                    titleLat = 0.1 * diff(latMap) + max(latMap);
                    switch ShowMapTitle, case true
                        m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' seasons.season{i} ')'], ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                    end
                    hold(mapAx, 'off')
                end
            end
            clearvars counter Dat samples nsamples lons lats labels lonMap latMap mapAx titleLon titleLat i j k l *ID
        end

        % Plots show between-cruise variability in sample locations.
        % It appears that some locations are sampled in multiple cruises; these are
        % probably the permanent mooring stations (Western Core Box?). Some sample
        % locations, however, appear only in a single cruise; these may cluster
        % around a permanent mooring station or form a transect. Two cruises
        % sampled a NW-SE transect southeast of South Georgia; the other cruises
        % sampled waters north & west of South Georgia, with some samples on the
        % north coast.

        % These data are non-trivial to organise!
        % First, group into 2 data sets: the southeast transect, and the
        % northwestern waters.
        % Then find coordinates of permanent mooring stations and group samples as
        % permanent stations or one-off casts.

        % Group cruises as 'northwest' or 'southeast' using sample longitudes
        % location = cell(cruises.ncruises, 1);
        MetaData.Location = cell(height(MetaData), 1);
        MetaData = movevars(MetaData, 'Location', 'after', 'Season');
        for i = 1:seasons.nseasons
            seasonID = seasons.fieldName{i};
            for j = 1:cruises.(seasonID).ncruises
                cruiseID = cruises.(seasonID).ID{j};
                ind = strcmp(MetaData.Season, seasons.season{i}) & ...
                    strcmp(MetaData.Cruise, cruiseID);
                if all(MetaData.Longitude(ind) < -36), r = 'northwest'; else, r = 'southeast'; end
                cruises.(seasonID).Location{j} = r;
                MetaData.Location(ind) = {r};
            end
        end
        clearvars ind r i j *ID

        % Aggregate data into these regional groups -- and convert Data struct into
        % a table, which should be easier to work with.
        for i = 1:seasons.nseasons
            if i == 1, clearvars Data2 cruises2; end
            seasonID = seasons.fieldName{i};
            for j = 1:cruises.(seasonID).ncruises
                cruiseID = cruises.(seasonID).ID{j};
                Dat = Data.(seasonID).(cruiseID);
                samples = fieldnames(Dat);
                for k = 1:length(samples)
                    sampleID = samples{k};
                    dat = Dat.(sampleID);
                    mdat = MetaData(strcmp(MetaData.Label, sampleID) & ...
                        strcmp(MetaData.Cruise, cruiseID),:);
                    fields = fieldnames(dat);
                    remove = fieldnames(mdat);
                    remove = remove(contains(remove, fields));
                    dat = rmfield(dat, remove);
                    mdat = repmat(mdat, [length(dat.Depth), 1]);
                    dat = struct2table(dat);
                    dat = [mdat, dat];
                    if i ~= 1 || j ~= 1 || k ~= 1
                        Data2 = [Data2; dat];
                    else
                        Data2 = dat;
                    end
                end
                dc = cruises.(seasonID);
                dc.ncruises = repmat(dc.ncruises, [1, dc.ncruises]);
                dc = struct2table(structfun(@(z) z(j), dc, 'UniformOutput', false));
                dc.Season = seasons.season(i);
                dc = movevars(dc, 'Season', 'before', 'ncruises');
                if i ~= 1 || j ~= 1
                    cruises2 = [cruises2; dc];
                else
                    cruises2 = dc;
                end
                %         if i == seasons.nseasons && j == cruises.(seasonID).ncruises
                %             cruises2.ncruises = repmat(cruises.ncruises, [height(cruises2), 1]);
                %         end
            end
        end
        Data = Data2;
        cruises = cruises2;
        clearvars Data2 cruises2 Dat dat mdat dc fields remove samples i j k *ID


        % The southeast data are a transect. The northwest data appear more
        % haphazard.
        % Organise the northwest data first...

        % Identify long term mooring sample stations and choose some labelling
        % scheme... Let's work from east to west, picking out obvious station groups
        Data.StationID = cell(height(Data), 1);
        Data = movevars(Data, 'StationID', 'after', 'Station');

        ndat = height(Data); % total number of measurements
        northwest = strcmp(Data.Location, 'northwest'); % index all data sampled from northwest
        southeast = strcmp(Data.Location, 'southeast'); % and southeast areas

        % 2 sample stations on South Georgia north coast (most distinctly plotted in 2014-2015, JR20141115)
        nstations = 2;
        whichCruise = strcmp(Data.Cruise, 'JR20141115');
        position = -55 < Data.Latitude & Data.Latitude < -54;
        ind = northwest & whichCruise & position; % indexes 2 sample stations from JR20141115
        % Find coordinates of these stations...
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind)); % index the first station
        ind_(:,2) = ind & ~ind_(:,1); % and second station
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % Distance (km) between these stations
        Dist = m_lldist(coords_(1,:), coords_(2,:));
        radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 1 and 2
        stationNames = {'1', '2'};
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether these stations were sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the 2 selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        % Cruise JR20120327 was unusual -- different pattern of sample stations
        indc = strcmp(cruises.ID, 'JR20120327'); % omit this cruise for now
        for m = 1:height(cruises)
            if indc(m), continue; end
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end
        % Now consider cruise JR20120327, where there's a dense cluster of samples
        % around station 1. Here, let's identify station 1, then assign other
        % samples names 1.1, 1.2, ...
        cruiseID = cruises.ID{indc};
        ind = position & strcmp(Data.Cruise, cruiseID);
        d = allDist(ind,:); % distance between indexed cruise samples and the 2 identified stations
        whichStation = any(d == min(d(:)));
        d = d(:,whichStation);
        % stationNames(whichStation)
        indm = ind & allDist(:,whichStation) == min(d);
        Data.StationID(indm) = stationNames(whichStation); % station 1
        ind_ = ind & cellfun(@(z) isempty(z), Data.StationID); % index all samples without StationID
        samples = unique(Data.Label(ind_), 'stable');
        % order samples from south to north
        d = unique(Data(ind_,{'Label', 'Longitude', 'Latitude'}));
        [~, o] = sort(d.Longitude, 'descend');
        d = d(o,:);
        [~, o] = sort(d.Latitude);
        d = d(o,:);
        for m = 1:length(samples)
            dm = d(m,:);
            indm = ind_ & strcmp(Data.Label, dm.Label);
            Data.StationID(indm) = {[stationNames{whichStation} '.' num2str(m)]};
        end


        % Now all samples from the northwest should have stations 1 and 2
        % identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % sum(strcmp(x, '1'))
        % sum(strcmp(x, '2'))

        % A second obvious cluster of 6 stations is off the northwest coast, most
        % distinctly plotted in cruises JR20141115 and JR15002 -- let's continue
        % using JR20141115 as the guide...
        nstations = 6;
        whichCruise = strcmp(Data.Cruise, 'JR20141115');
        position = -54 < Data.Latitude & Data.Latitude < -53 & ...
            -42 < Data.Longitude & Data.Longitude < -36;
        ind = northwest & whichCruise & position; % indexes 6 sample stations from JR20141115
        % Find coordinates of these stations... ordered east to west then south to north
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,2) = ind & Data.Longitude < max(Data.Longitude(ind)) & Data.Latitude < unique(Data.Latitude(ind_(:,1)));
        ind_(:,3) = ind_(:,2) & Data.Longitude < max(Data.Longitude(ind_(:,2)));
        ind_(:,2) = ind_(:,2) & Data.Longitude > min(Data.Longitude(ind_(:,2)));
        ind_(:,4) = ind & Data.Latitude > max(Data.Latitude(ind_(:,1)));
        ind_(:,6) = ind_(:,4) & Data.Longitude == max(Data.Longitude(ind_(:,4)));
        ind_(:,5) = ind_(:,4) & ~ind_(:,6);
        ind_(:,4) = ind_(:,4) & Data.Longitude == min(Data.Longitude(ind_(:,4)));
        ind_(:,5) = ind_(:,5) & Data.Longitude == max(Data.Longitude(ind_(:,5)));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 3-8
        stationNames = cellstr(string(3:8));
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether these stations were sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the 6 selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Now all samples from the northwest should have stations 3-8 identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % [sum(strcmp(x, '3')) , sum(strcmp(x, '4')), sum(strcmp(x, '5')) , sum(strcmp(x, '6')), sum(strcmp(x, '7')) , sum(strcmp(x, '8'))]

        % Once again, cruise JR20120327 is unusual as the sample does not exactly
        % correspond to a station location.
        % Find the nearest station (x) to the JR20131112 samples then name this x.1
        cruiseID = 'JR20120327';
        ind = strcmp(Data.Cruise, cruiseID);
        d = allDist(ind,:); % distance between all cruise samples and the 6 identified stations
        whichStation = any(d == min(d(:))); % nearest station
        d = d(:,whichStation);
        indm = ind & allDist(:,whichStation) == min(d);
        Data.StationID(indm) = {[stationNames{whichStation} '.1']};



        % The third obvious cluster of 2 stations is further northwest, most
        % distinctly plotted in cruises JR15002.
        nstations = 2;
        whichCruise = strcmp(Data.Cruise, 'JR15002');
        position = Data.Latitude > -53 & Data.Longitude > -41;
        ind = northwest & whichCruise & position; % indexes 2 sample stations from JR15002
        % Find coordinates of these stations... ordered east to west
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,2) = ind & Data.Longitude == min(Data.Longitude(ind));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 9 and 10
        stationNames = cellstr(string(9:10));
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether these stations were sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the 2 selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Now all samples from the northwest should have stations 9 and 10 identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % [sum(strcmp(x, '9')) , sum(strcmp(x, '10'))]

        % One sample from JR20141115 does not correspond to a station, but is
        % nearby station 10 so label it as station 10.1.
        cruiseID = 'JR20141115';
        ind = strcmp(Data.Cruise, cruiseID) & position;
        ind_ = ind & cellfun(@(z) isempty(z), Data.StationID);
        d = allDist(ind_,:);
        whichStation = any(d == min(d(:)));
        Data.StationID(ind_) = {[stationNames{whichStation}, '.1']};


        % The fourth (and final) station group from the northwest samples is
        % southwest of South Georgia -- use cruise JR20141115 as the guide
        nstations = 1;
        whichCruise = strcmp(Data.Cruise, 'JR20141115');
        position = Data.Latitude < -54 & Data.Longitude < -39;
        ind = northwest & whichCruise & position; % indexes sample station from JR20141115
        % Find coordinates of this station
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define this as station 11
        stationNames = {'11'};
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether this station was sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the selected station
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Now all samples from the northwest should have stations 9 and 10 identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % [sum(strcmp(x, '11'))]

        % Cruise JR16003 contains 3 extra samples forming an eastward transect from
        % the identified station 11 -- label these as 11.1-11.3.
        cruiseID = 'JR16003';
        ind = strcmp(Data.Cruise, cruiseID) & position;
        ind = ind & cellfun(@(z) isempty(z), Data.StationID); % index all samples to label
        ind_ = false(ndat, 3);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,3) = ind & Data.Longitude == min(Data.Longitude(ind));
        ind_(:,2) = ind & ~(ind_(:,1) | ind_(:,3));
        for m = 1:3, Data.StationID(ind_(:,m)) = {[stationNames{:}, '.' num2str(m)]}; end


        % Cruise JR20120327 has 2 extra samples not associated to other stations --
        % label these as stations 12 and 13
        nstations = 2;
        whichCruise = strcmp(Data.Cruise, 'JR20120327');
        position = Data.Longitude < -41 & Data.Latitude > -54;
        ind = northwest & whichCruise & position; % indexes sample station from JR20141115
        % Find coordinates of these stations
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,2) = ind & Data.Longitude == min(Data.Longitude(ind));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 12 and 13
        stationNames = cellstr(string(12:13));
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether this station was sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the selected station
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Update the MetaData
        MetaData.StationID = cell(height(MetaData), 1);
        MetaData = movevars(MetaData, 'StationID', 'after', 'Station');
        for i = 1:height(MetaData)
            if ~strcmp(MetaData.Location{i}, 'northwest'), continue; end
            ind = strcmp(Data.Cruise, MetaData.Cruise{i}) & strcmp(Data.Label, MetaData.Label{i});
            MetaData.StationID(i) = unique(Data.StationID(ind));
        end

        % Replot the maps for the northwest cruises
        if plotMaps
            close(plotHandles.map2); close(plotHandles.map3); close(plotHandles.map4); close(plotHandles.map6); close(plotHandles.map7);
            for i = 1:height(cruises)
                if ~strcmp(cruises.Location{i}, 'northwest'), continue; end
                cruiseID = cruises.ID{i};
                md = MetaData(strcmp(MetaData.Cruise, cruiseID),:);
                lons = md.Longitude;
                lats = md.Latitude;
                labels = md.StationID;
                % Create South Georgia base map
                mapName = ['map' num2str(i)];
                figure
                plotHandles.(mapName) = gcf;
                %         assignin('base', mapName, figure)
                plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                    'XaxisLocation', XaxisLocation);
                % Extract map bounding coordinates
                [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0]);
                mapAx = gca;
                % Plot sample positions
                hold(mapAx, 'on')
                m_scatter(lons, lats, 18, [1 0 0], 'filled')
                % remove duplicate labels -- multiple samples at single station
                lt = table(lons, lats, labels);
                nlt = cell2mat(cellfun(@(z) str2num(z), lt.labels, 'UniformOutput', false));
                [~, o] = sort(nlt);
                lt = lt(o,:); nlt = nlt(o,:);
                ult = unique(nlt);
                for j = 1:length(ult)
                    k = nlt == ult(j);
                    x = lt(find(k,1),:);
                    lt(k,:) = []; nlt(k) = [];
                    lt = [lt; x]; nlt = [nlt; ult(j)];
                end
                m_text(lt.lons, lt.lats, lt.labels, 'FontSize', 8, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                titleLon = mean(lonMap);
                titleLat = 0.1 * diff(latMap) + max(latMap);
                switch ShowMapTitle, case true
                    m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' cruises.Season{i} ')'], ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                end
                hold(mapAx, 'off')
            end
        end

        % Station labels are now consistent between cruises sampling the northwest
        % region... now organise the southeast region

        % Use the 9 samples from cruise JR20150309 to set the station IDs.
        nstations = 9;
        stationNames = cellstr(string(1:nstations));

        cruiseID = 'JR20150309';
        whichCruise = strcmp(Data.Cruise, cruiseID);
        d = Data(whichCruise,:);
        samples = unique(d.Label);
        for i = 1:length(samples)
            sampleID = samples{i};
            whichSample = strcmp(d.Label, sampleID);
            d.StationID(whichSample) = repmat(stationNames(i), [sum(whichSample), 1]);
        end
        Data(whichCruise,:) = d;

        % position = Data.Longitude < -41 & Data.Latitude > -54;
        ind = southeast & whichCruise;
        % Find coordinates of these stations -- listed south -> north
        ind_ = false(ndat, nstations);
        for i = 1:nstations, ind_(:,i) = ind & strcmp(Data.StationID, num2str(i)); end
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % Get station IDs for other southeast cruises
        % Find distances between all samples and the selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Update the MetaData
        for i = 1:height(MetaData)
            if ~strcmp(MetaData.Location{i}, 'southeast'), continue; end
            ind = strcmp(Data.Cruise, MetaData.Cruise{i}) & strcmp(Data.Label, MetaData.Label{i});
            MetaData.StationID(i) = unique(Data.StationID(ind));
        end

        % Replot the maps for the southeast cruises
        if plotMaps
            %     delete(map1); delete(map5)
            %     close('map1', 'map5')
            %     close(map1); close(map5)
            close(plotHandles.map1); close(plotHandles.map5);
            for i = 1:height(cruises)
                if ~strcmp(cruises.Location{i}, 'southeast'), continue; end
                cruiseID = cruises.ID{i};
                md = MetaData(strcmp(MetaData.Cruise, cruiseID),:);
                lons = md.Longitude;
                lats = md.Latitude;
                labels = md.StationID;
                % Create South Georgia base map
                mapName = ['map' num2str(i)];
                figure
                plotHandles.(mapName) = gcf;
                %         assignin('base', mapName, figure)
                plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                    'XaxisLocation', XaxisLocation);
                % Extract map bounding coordinates
                [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0]);
                mapAx = gca;
                % Plot sample positions
                hold(mapAx, 'on')
                m_scatter(lons, lats, 18, [1 0 0], 'filled')
                % remove duplicate labels -- multiple samples at single station
                lt = table(lons, lats, labels);
                nlt = cell2mat(cellfun(@(z) str2double(z), lt.labels, 'UniformOutput', false));
                [~, o] = sort(nlt);
                lt = lt(o,:); nlt = nlt(o,:);
                ult = unique(nlt);
                for j = 1:length(ult)
                    k = nlt == ult(j);
                    x = lt(find(k,1),:);
                    lt(k,:) = []; nlt(k) = [];
                    lt = [lt; x]; nlt = [nlt; ult(j)];
                end
                m_text(lt.lons, lt.lats, lt.labels, 'FontSize', 8, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                titleLon = mean(lonMap);
                titleLat = 0.1 * diff(latMap) + max(latMap);
                switch ShowMapTitle, case true
                    m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' cruises.Season{i} ')'], ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                end
                hold(mapAx, 'off')
            end
        end


        %% Are samples depth profiles?

        % Not all of the samples are depth profiles -- many appear to be
        % measurements at/around a single depth, often near the seafloor.
        % Determine which sample are depth profiles and index them.

        % This can be automated using coefficient of variation of measurement depth
        % -- very small CoV => measurements cluster around single depth.

        % CoV_lim = 0.05; % CoV >= CoV_lim => depth profile, CoV < CoV_lim => sample around single depth. CoV_lim is arbitrary choice...
        Data.DepthCoV = nan(height(Data),1);
        Data.DepthProfile = false(height(Data),1);
        MetaData.DepthCoV = false(height(MetaData), 1);
        MetaData.DepthProfile = false(height(MetaData), 1);

        for i = 1:height(cruises)
            whichCruise = cruises.ID{i};
            indi = strcmp(Data.Cruise, whichCruise); % index cruise
            allStations = unique(Data.StationID(indi), 'stable');
            for j = 1:length(allStations)
                whichStation = allStations{j};
                indj = indi & strcmp(Data.StationID, whichStation);  % index cruise & station
                allSamples = unique(Data.Label(indj));
                for k = 1:length(allSamples)
                    whichSample = allSamples{k};
                    indk = indj & strcmp(Data.Label, whichSample);
                    indm = strcmp(MetaData.Cruise, whichCruise) & strcmp(MetaData.StationID, whichStation) & strcmp(MetaData.Label, whichSample);
                    CoV = std(Data.Depth(indk), 'omitnan') / mean(Data.Depth(indk), 'omitnan');
                    Data.DepthCoV(indk) = CoV;
                    MetaData.DepthCoV(indm) = CoV;
                    if CoV >= CoV_lim
                        Data.DepthProfile(indk) = true;
                        MetaData.DepthProfile(indm) = true;
                    end
                    %             if i ~= 1 || j ~= 1 || k ~= 1
                    %                 profileSummary.Cruise = [profileSummary.Cruise; whichCruise];
                    %                 profileSummary.StationID = [profileSummary.StationID; whichStation];
                    %                 profileSummary.Label = [profileSummary.Label; whichSample];
                    %                 profileSummary.DepthCoV = [profileSummary.DepthCoV; CoV];
                    %                 profileSummary.DepthProfile = [profileSummary.DepthProfile; CoV >= CoV_lim];
                    %             else
                    %                 profileSummary.Cruise = {whichCruise};
                    %                 profileSummary.StationID = {whichStation};
                    %                 profileSummary.Label = {whichSample};
                    %                 profileSummary.DepthCoV = CoV;
                    %                 profileSummary.DepthProfile = CoV >= CoV_lim;
                    %             end
                end
            end
        end

        % profileSummary = struct2table(profileSummary);

        t = sprintf('%u%%', round(100 * sum(MetaData.DepthProfile) / height(MetaData)));
        t = sprintf('Of all CTD samples, %s', t);
        fprintf('\n\n'); fprintf(1, '%s are depth profiles.\nThe remaining samples are measurements around a single depth.\n(Automatic determination of depth profiles depends on CoV_lim,\nwhich is a somewhat arbitrary choice so maybe try ajusting\nthis parameter... lower values include more samples.)', t); fprintf('\n\n')


        %% Flag spurious points

        % Although the quality control flags suggest all measurements are valid,
        % plotting demonstrates lots of spurious measurements that mostly appear
        % in the first few data points. This suggests that sensor readings are
        % unreliable when CTD is first deployed, so these measurements should be
        % omitted.
        % Focus on depth profile samples -- ignoring all others for now.

        % Each sample is a series of successive measurements seprated by some,
        % hopefully uniform, time increment.
        % Plot measurement depth against measurement number -- separate multipanel
        % plot for each cruise.

        if plotSpuriousPoints
            for i = 1:height(cruises)
                whichCruise = cruises.ID{i};
                ps = MetaData(strcmp(MetaData.Cruise, whichCruise),:);
                if ~any(ps.DepthProfile), continue; end
                np = sum(ps.DepthProfile);
                nrows = ceil(np ^ 0.5); ncols = ceil(np / nrows);
                testPlotHandles.(['plot' num2str(i)]) = figure;
                set(gcf, {'Unit','Position'}, {'inches', [0 0 4*ncols 3*nrows]})
                j = 0;
                for k = 1:height(ps)
                    if ~ps.DepthProfile(k), continue; end
                    j = j + 1;
                    subplot(nrows, ncols, j)
                    whichStation = ps.StationID{k};
                    whichSample = ps.Label{k};
                    ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation) & strcmp(Data.Label, whichSample);
                    Dat = Data(ind,:);
                    sampleNumber = 1:height(Dat);
                    scatter(sampleNumber, -Dat.Depth)
                    xlabel('Measurement'); ylabel('Depth (m)'); title([whichCruise ' - station ' whichStation ' - ' whichSample], 'FontWeight', 'normal')
                end
            end
        end

        % Clearly, many depth profiles contain spurious values, mostly appearing as
        % the initial measurements.
        % Points forming lines with negative gradient correspond to good
        % measurements along a depth profile.

        % The following is an ad-hoc method to catch spurious points... it's imperfect

        Data.SpuriousMeasure = false(height(Data),1);

        % This is a decent method, but it doesn't quite catch all the spurious
        % points... Just manually find the spurious points, using the following
        % plots as a guide...

        for i = 1:height(cruises)
            whichCruise = cruises.ID{i};
            indi = strcmp(Data.Cruise, whichCruise); % index cruise
            allStations = cellstr(string(unique(cellfun(@(z) str2double(z), Data.StationID(indi)))));
            for j = 1:length(allStations)
                whichStation = allStations{j};
                indj = indi & strcmp(Data.StationID, whichStation); % index cruise & station
                allSamples = unique(Data.Label(indj));
                for k = 1:length(allSamples)
                    whichSample = allSamples{k};
                    indk = indj & strcmp(Data.Label, whichSample); % index cruise & station & sample
                    isDepthProfile = MetaData.DepthProfile(strcmp(MetaData.Cruise, whichCruise) & strcmp(MetaData.StationID, whichStation) & strcmp(MetaData.Label, whichSample));
                    if ~isDepthProfile, continue; end
                    Dat = Data(indk,:);
                    np = height(Dat); % number of measurements
                    mp = 1:np;
                    badPoints = false(np,1); % index spurious data points
                    %             figure
                    %             scatter(1:sum(indk), -Dat.Depth); xlabel('Measurement'); ylabel('Depth (m)')

                    % Are the final measurements OK? The last point should be the
                    % deepest. If it's not then there's a problem...
                    deepestPoint = Dat.Depth == max(Dat.Depth);
                    if ~deepestPoint(end)
                        badPoints = badPoints | logical(cumsum([0; deepestPoint(1:end-1)]));
                    end
                    nbp = sum(badPoints);

                    %             hold on
                    %             scatter(mp(badPoints), -Dat.Depth(badPoints), 'MarkerEdgeColor', [1 0 0])
                    %             hold off
                    % Spurious final measurements have (hopefully) been identified,
                    % but note that this method will fail if the deepest samples
                    % appear as spurious initial measurements... it looks OK for
                    % the BODC South Georgia.
                    % Now remove spurious initial measurements.
                    d = flip(Dat.Depth);
                    badPoints = flip(badPoints);
                    %             scatter(mp, d) % isolate points on the linear slope of negative gradient
                    %             hold on
                    %             scatter(mp(badPoints), d(badPoints), 'MarkerEdgeColor', [1 0 0])
                    %             hold off
                    d_ = d(~badPoints); % omit any spurious final points
                    dd = diff(d_); dd = [dd(1); dd];
                    ma = ceil(movstdRange * np);
                    mm = movmean(dd, [ma, 0]);
                    ms = movstd(dd, [ma, 0]);
                    mm = [mm(1); mm(1:end-1)]; % shift by one to get expected values
                    ms = [ms(1); ms(1:end-1)];
                    %             ms(1:nbp+1) = ms(nbp+2); % remove zeros

                    cutOffMetric = zeros(np, leeway+1);
                    cutOffMetric(nbp+1:end,1) = abs(ms ./ dd); % .^ 2;
                    for l = 1:leeway
                        cutOffMetric(nbp+1:end,l+1) = abs([ms(1:end-l) ./ dd(l+1:end); zeros(l, 1)]);
                    end
                    cutOffMetric(1:ma,:) = 0; % avoid false positives at final measurements
                    badPoints = badPoints | cumsum(all(cutOffMetric > cutOffValue, 2)) > 0;
                    badPoints = flip(badPoints);
                    %             scatter(mp, -Dat.Depth); xlabel('Measurement'); ylabel('Depth (m)')
                    %             hold on
                    %             scatter(mp(badPoints), -Dat.Depth(badPoints), 'MarkerEdgeColor', [1 0 0])
                    %             hold off
                    Data.SpuriousMeasure(indk) = badPoints;
                end
            end
        end

        if plotSpuriousPoints
            % Replot to view spurious measurements
            fields = fieldnames(testPlotHandles);
            for i = 1:length(fields), close(testPlotHandles.(fields{i})); end
            clearvars testPlotHandles

            for i = 1:height(cruises)
                whichCruise = cruises.ID{i};
                ps = MetaData(strcmp(MetaData.Cruise, whichCruise),:);
                if ~any(ps.DepthProfile), continue; end
                np = sum(ps.DepthProfile);
                nrows = ceil(np ^ 0.5); ncols = ceil(np / nrows);
                testPlotHandles.(['plot' num2str(i)]) = figure;
                set(gcf, {'Unit','Position'}, {'inches', [0 0 4*ncols 3*nrows]})
                j = 0;
                for k = 1:height(ps)
                    if ~ps.DepthProfile(k), continue; end
                    j = j + 1;
                    subplot(nrows, ncols, j)
                    whichStation = ps.StationID{k};
                    whichSample = ps.Label{k};
                    ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation) & strcmp(Data.Label, whichSample);
                    Dat = Data(ind,:);
                    sampleNumber = 1:height(Dat);
                    scatter(sampleNumber, -Dat.Depth)
                    xlabel('Measurement'); ylabel('Depth (m)'); title([whichCruise ' - station ' whichStation ' - ' whichSample], 'FontWeight', 'normal')
                    hold on
                    spuriousPoints = Dat.SpuriousMeasure;
                    scatter(sampleNumber(spuriousPoints), -Dat.Depth(spuriousPoints), 'MarkerEdgeColor', [1 0 0])
                    hold off
                end
            end
        end

        % Now, with help from these plots, manually specify the spurious points...
        switch source
            case 'BODC'
                whichCruise = 'JR20120327'; % no spurious points
                Data.SpuriousMeasure(strcmp(Data.Cruise, whichCruise)) = false;
                whichCruise = 'JR20131112'; % the above filtering method worked OK
                whichCruise = 'JR20141115'; % stations 3 & 4 need fixed
                whichStation = '3';
                ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation);
                x = false(sum(ind), 1); x(1:16) = true;
                Data.SpuriousMeasure(ind) = x;
                whichStation = '4';
                ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation);
                x = false(sum(ind), 1); x(1:57) = true;
                Data.SpuriousMeasure(ind) = x;
                whichCruise = 'JR15002'; % the above filtering method worked OK
                whichCruise = 'JR16003'; % the above filtering method worked OK
        end


        if plotSpuriousPoints

            % Replot to view spurious measurements
            fields = fieldnames(testPlotHandles);
            for i = 1:length(fields), close(testPlotHandles.(fields{i})); end
            clearvars testPlotHandles

            for i = 1:height(cruises)
                whichCruise = cruises.ID{i};
                ps = MetaData(strcmp(MetaData.Cruise, whichCruise),:);
                if ~any(ps.DepthProfile), continue; end
                np = sum(ps.DepthProfile);
                nrows = ceil(np ^ 0.5); ncols = ceil(np / nrows);
                testPlotHandles.(['plot' num2str(i)]) = figure;
                set(gcf, {'Unit','Position'}, {'inches', [0 0 4*ncols 3*nrows]})
                j = 0;
                for k = 1:height(ps)
                    if ~ps.DepthProfile(k), continue; end
                    j = j + 1;
                    subplot(nrows, ncols, j)
                    whichStation = ps.StationID{k};
                    whichSample = ps.Label{k};
                    ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation) & strcmp(Data.Label, whichSample);
                    Dat = Data(ind,:);
                    sampleNumber = 1:height(Dat);
                    scatter(sampleNumber, -Dat.Depth)
                    xlabel('Measurement'); ylabel('Depth (m)'); title([whichCruise ' - station ' whichStation ' - ' whichSample], 'FontWeight', 'normal')
                    hold on
                    spuriousPoints = Dat.SpuriousMeasure;
                    scatter(sampleNumber(spuriousPoints), -Dat.Depth(spuriousPoints), 'MarkerEdgeColor', [1 0 0])
                    hold off
                end
                sgtitle('Depth profiles with spurious points highlighted')
            end

        end


        %% Tidy up

        % Check for and omit any negative (above surface) depths
        Data = Data(Data.Depth > 0,:);

        % Remove CTD casts that are not depth profiles?
        switch onlyDepthProfiles, case true
            Data = Data(Data.DepthProfile,:);
            MetaData = MetaData(MetaData.DepthProfile,:);
            fprintf('\n\n'); fprintf('The data have been filtered to only include CTD depth profiles.\nAll other CTD casts are excluded and not shown on maps.'); fprintf('\n\n')

            % Replot maps to show only depth profiles
            if plotMaps
                fields = fieldnames(plotHandles);
                for i = 1:length(fields), close(plotHandles.(fields{i})); end
                clearvars plotHandles

                % Create a separate map for each cruise
                ucruises = unique(MetaData.Cruise, 'stable');
                ncruises = length(ucruises);
                for i = 1:ncruises
                    cruiseID = ucruises{i};
                    seasonID = cruises.Season{strcmp(cruises.ID, cruiseID)};
                    Dat = Data(strcmp(Data.Cruise, cruiseID),:);
                    mDat = MetaData(strcmp(MetaData.Cruise, cruiseID),:);
                    labels = mDat.StationID;
                    %             samples = mDat.Label;
                    nsamples = length(labels);
                    lons = mDat.Longitude; lats = mDat.Latitude;
                    %             labels = cellfun(@(z) num2str(str2double(z(4:end))), mDat.Label, 'UniformOutput', false);

                    % Create South Georgia base map
                    mapName = ['map' num2str(i)];
                    figure
                    plotHandles.(mapName) = gcf;
                    %             assignin('base', mapName, figure)
                    plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                        'XaxisLocation', XaxisLocation);
                    % Extract map bounding coordinates
                    [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0]);
                    mapAx = gca;
                    % Plot sample positions
                    hold(mapAx, 'on')
                    m_scatter(lons, lats, 18, [1 0 0], 'filled')
                    m_text(lons, lats, labels, 'FontSize', 8, ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                    titleLon = mean(lonMap);
                    titleLat = 0.1 * diff(latMap) + max(latMap);
                    switch ShowMapTitle, case true
                        m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' seasonID ')'], ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                    end
                    hold(mapAx, 'off')
                end
                clearvars counter Dat samples nsamples lons lats labels lonMap latMap mapAx titleLon titleLat i j k l *ID
            end
        end


        % Remove spurious points?
        switch removeSpuriousPoints
            case true
                Data = Data(~Data.SpuriousMeasure,:);
            case false
                disp('Spurious measurements may remain in the data. Examine "Data.SpuriousMeasure" to check.')
        end



    case 'BAS central storage'
        %% These data should already have been cleaned, so just make maps
        coordsTable = readtable('regional bounding coordinates.csv');

        if plotMaps
            figure
            plotBaseMap('Southern Ocean', 'createMap', true, ...
                'coordsTable', coordsTable, 'redrawCoastline', true, ...
                'edgecolour', [0 0 0]);

            Data
        end





        if plotMaps
            % Create a separate map for each cruise
            for i = 1:seasons.nseasons
                if i == 1, counter = 0; end
                seasonID = seasons.fieldName{i};
                for j = 1:cruises.(seasonID).ncruises
                    counter = counter + 1;
                    cruiseID = cruises.(seasonID).ID{j};
                    Dat = Data.(seasonID).(cruiseID);
                    samples = fieldnames(Dat);
                    nsamples = length(samples);
                    % Extract CTD sample locations & create ID tags
                    lons = nan(nsamples, 1); lats = nan(nsamples, 1); labels = cell(nsamples, 1);
                    for k = 1:nsamples
                        sampleID = samples{k};
                        lons(k) = Dat.(sampleID).Longitude;
                        lats(k) = Dat.(sampleID).Latitude;
                        l = Dat.(sampleID).Label;
                        labels{k} = num2str(str2double(l(4:end)));
                    end
                    % Create South Georgia base map
                    mapName = ['map' num2str(counter)];
                    figure
                    plotHandles.(mapName) = gcf;
                    %             assignin('base', mapName, figure)
                    plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                        'XaxisLocation', XaxisLocation);
                    % Extract map bounding coordinates
                    [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0]);
                    mapAx = gca;
                    % Plot sample positions
                    hold(mapAx, 'on')
                    m_scatter(lons, lats, 18, [1 0 0], 'filled')
                    m_text(lons, lats, labels, 'FontSize', 8, ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                    titleLon = mean(lonMap);
                    titleLat = 0.1 * diff(latMap) + max(latMap);
                    switch ShowMapTitle, case true
                        m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' seasons.season{i} ')'], ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                    end
                    hold(mapAx, 'off')
                end
            end
            clearvars counter Dat samples nsamples lons lats labels lonMap latMap mapAx titleLon titleLat i j k l *ID
        end

        % Plots show between-cruise variability in sample locations.
        % It appears that some locations are sampled in multiple cruises; these are
        % probably the permanent mooring stations (Western Core Box?). Some sample
        % locations, however, appear only in a single cruise; these may cluster
        % around a permanent mooring station or form a transect. Two cruises
        % sampled a NW-SE transect southeast of South Georgia; the other cruises
        % sampled waters north & west of South Georgia, with some samples on the
        % north coast.

        % These data are non-trivial to organise!
        % First, group into 2 data sets: the southeast transect, and the
        % northwestern waters.
        % Then find coordinates of permanent mooring stations and group samples as
        % permanent stations or one-off casts.

        % Group cruises as 'northwest' or 'southeast' using sample longitudes
        % location = cell(cruises.ncruises, 1);
        MetaData.Location = cell(height(MetaData), 1);
        MetaData = movevars(MetaData, 'Location', 'after', 'Season');
        for i = 1:seasons.nseasons
            seasonID = seasons.fieldName{i};
            for j = 1:cruises.(seasonID).ncruises
                cruiseID = cruises.(seasonID).ID{j};
                ind = strcmp(MetaData.Season, seasons.season{i}) & ...
                    strcmp(MetaData.Cruise, cruiseID);
                if all(MetaData.Longitude(ind) < -36), r = 'northwest'; else, r = 'southeast'; end
                cruises.(seasonID).Location{j} = r;
                MetaData.Location(ind) = {r};
            end
        end
        clearvars ind r i j *ID

        % Aggregate data into these regional groups -- and convert Data struct into
        % a table, which should be easier to work with.
        for i = 1:seasons.nseasons
            if i == 1, clearvars Data2 cruises2; end
            seasonID = seasons.fieldName{i};
            for j = 1:cruises.(seasonID).ncruises
                cruiseID = cruises.(seasonID).ID{j};
                Dat = Data.(seasonID).(cruiseID);
                samples = fieldnames(Dat);
                for k = 1:length(samples)
                    sampleID = samples{k};
                    dat = Dat.(sampleID);
                    mdat = MetaData(strcmp(MetaData.Label, sampleID) & ...
                        strcmp(MetaData.Cruise, cruiseID),:);
                    fields = fieldnames(dat);
                    remove = fieldnames(mdat);
                    remove = remove(contains(remove, fields));
                    dat = rmfield(dat, remove);
                    mdat = repmat(mdat, [length(dat.Depth), 1]);
                    dat = struct2table(dat);
                    dat = [mdat, dat];
                    if i ~= 1 || j ~= 1 || k ~= 1
                        Data2 = [Data2; dat];
                    else
                        Data2 = dat;
                    end
                end
                dc = cruises.(seasonID);
                dc.ncruises = repmat(dc.ncruises, [1, dc.ncruises]);
                dc = struct2table(structfun(@(z) z(j), dc, 'UniformOutput', false));
                dc.Season = seasons.season(i);
                dc = movevars(dc, 'Season', 'before', 'ncruises');
                if i ~= 1 || j ~= 1
                    cruises2 = [cruises2; dc];
                else
                    cruises2 = dc;
                end
                %         if i == seasons.nseasons && j == cruises.(seasonID).ncruises
                %             cruises2.ncruises = repmat(cruises.ncruises, [height(cruises2), 1]);
                %         end
            end
        end
        Data = Data2;
        cruises = cruises2;
        clearvars Data2 cruises2 Dat dat mdat dc fields remove samples i j k *ID


        % The southeast data are a transect. The northwest data appear more
        % haphazard.
        % Organise the northwest data first...

        % Identify long term mooring sample stations and choose some labelling
        % scheme... Let's work from east to west, picking out obvious station groups
        Data.StationID = cell(height(Data), 1);
        Data = movevars(Data, 'StationID', 'after', 'Station');

        ndat = height(Data); % total number of measurements
        northwest = strcmp(Data.Location, 'northwest'); % index all data sampled from northwest
        southeast = strcmp(Data.Location, 'southeast'); % and southeast areas

        % 2 sample stations on South Georgia north coast (most distinctly plotted in 2014-2015, JR20141115)
        nstations = 2;
        whichCruise = strcmp(Data.Cruise, 'JR20141115');
        position = -55 < Data.Latitude & Data.Latitude < -54;
        ind = northwest & whichCruise & position; % indexes 2 sample stations from JR20141115
        % Find coordinates of these stations...
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind)); % index the first station
        ind_(:,2) = ind & ~ind_(:,1); % and second station
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % Distance (km) between these stations
        Dist = m_lldist(coords_(1,:), coords_(2,:));
        radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 1 and 2
        stationNames = {'1', '2'};
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether these stations were sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the 2 selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        % Cruise JR20120327 was unusual -- different pattern of sample stations
        indc = strcmp(cruises.ID, 'JR20120327'); % omit this cruise for now
        for m = 1:height(cruises)
            if indc(m), continue; end
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end
        % Now consider cruise JR20120327, where there's a dense cluster of samples
        % around station 1. Here, let's identify station 1, then assign other
        % samples names 1.1, 1.2, ...
        cruiseID = cruises.ID{indc};
        ind = position & strcmp(Data.Cruise, cruiseID);
        d = allDist(ind,:); % distance between indexed cruise samples and the 2 identified stations
        whichStation = any(d == min(d(:)));
        d = d(:,whichStation);
        % stationNames(whichStation)
        indm = ind & allDist(:,whichStation) == min(d);
        Data.StationID(indm) = stationNames(whichStation); % station 1
        ind_ = ind & cellfun(@(z) isempty(z), Data.StationID); % index all samples without StationID
        samples = unique(Data.Label(ind_), 'stable');
        % order samples from south to north
        d = unique(Data(ind_,{'Label', 'Longitude', 'Latitude'}));
        [~, o] = sort(d.Longitude, 'descend');
        d = d(o,:);
        [~, o] = sort(d.Latitude);
        d = d(o,:);
        for m = 1:length(samples)
            dm = d(m,:);
            indm = ind_ & strcmp(Data.Label, dm.Label);
            Data.StationID(indm) = {[stationNames{whichStation} '.' num2str(m)]};
        end


        % Now all samples from the northwest should have stations 1 and 2
        % identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % sum(strcmp(x, '1'))
        % sum(strcmp(x, '2'))

        % A second obvious cluster of 6 stations is off the northwest coast, most
        % distinctly plotted in cruises JR20141115 and JR15002 -- let's continue
        % using JR20141115 as the guide...
        nstations = 6;
        whichCruise = strcmp(Data.Cruise, 'JR20141115');
        position = -54 < Data.Latitude & Data.Latitude < -53 & ...
            -42 < Data.Longitude & Data.Longitude < -36;
        ind = northwest & whichCruise & position; % indexes 6 sample stations from JR20141115
        % Find coordinates of these stations... ordered east to west then south to north
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,2) = ind & Data.Longitude < max(Data.Longitude(ind)) & Data.Latitude < unique(Data.Latitude(ind_(:,1)));
        ind_(:,3) = ind_(:,2) & Data.Longitude < max(Data.Longitude(ind_(:,2)));
        ind_(:,2) = ind_(:,2) & Data.Longitude > min(Data.Longitude(ind_(:,2)));
        ind_(:,4) = ind & Data.Latitude > max(Data.Latitude(ind_(:,1)));
        ind_(:,6) = ind_(:,4) & Data.Longitude == max(Data.Longitude(ind_(:,4)));
        ind_(:,5) = ind_(:,4) & ~ind_(:,6);
        ind_(:,4) = ind_(:,4) & Data.Longitude == min(Data.Longitude(ind_(:,4)));
        ind_(:,5) = ind_(:,5) & Data.Longitude == max(Data.Longitude(ind_(:,5)));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 3-8
        stationNames = cellstr(string(3:8));
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether these stations were sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the 6 selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Now all samples from the northwest should have stations 3-8 identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % [sum(strcmp(x, '3')) , sum(strcmp(x, '4')), sum(strcmp(x, '5')) , sum(strcmp(x, '6')), sum(strcmp(x, '7')) , sum(strcmp(x, '8'))]

        % Once again, cruise JR20120327 is unusual as the sample does not exactly
        % correspond to a station location.
        % Find the nearest station (x) to the JR20131112 samples then name this x.1
        cruiseID = 'JR20120327';
        ind = strcmp(Data.Cruise, cruiseID);
        d = allDist(ind,:); % distance between all cruise samples and the 6 identified stations
        whichStation = any(d == min(d(:))); % nearest station
        d = d(:,whichStation);
        indm = ind & allDist(:,whichStation) == min(d);
        Data.StationID(indm) = {[stationNames{whichStation} '.1']};



        % The third obvious cluster of 2 stations is further northwest, most
        % distinctly plotted in cruises JR15002.
        nstations = 2;
        whichCruise = strcmp(Data.Cruise, 'JR15002');
        position = Data.Latitude > -53 & Data.Longitude > -41;
        ind = northwest & whichCruise & position; % indexes 2 sample stations from JR15002
        % Find coordinates of these stations... ordered east to west
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,2) = ind & Data.Longitude == min(Data.Longitude(ind));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 9 and 10
        stationNames = cellstr(string(9:10));
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether these stations were sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the 2 selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Now all samples from the northwest should have stations 9 and 10 identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % [sum(strcmp(x, '9')) , sum(strcmp(x, '10'))]

        % One sample from JR20141115 does not correspond to a station, but is
        % nearby station 10 so label it as station 10.1.
        cruiseID = 'JR20141115';
        ind = strcmp(Data.Cruise, cruiseID) & position;
        ind_ = ind & cellfun(@(z) isempty(z), Data.StationID);
        d = allDist(ind_,:);
        whichStation = any(d == min(d(:)));
        Data.StationID(ind_) = {[stationNames{whichStation}, '.1']};


        % The fourth (and final) station group from the northwest samples is
        % southwest of South Georgia -- use cruise JR20141115 as the guide
        nstations = 1;
        whichCruise = strcmp(Data.Cruise, 'JR20141115');
        position = Data.Latitude < -54 & Data.Longitude < -39;
        ind = northwest & whichCruise & position; % indexes sample station from JR20141115
        % Find coordinates of this station
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define this as station 11
        stationNames = {'11'};
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether this station was sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the selected station
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Now all samples from the northwest should have stations 9 and 10 identified...
        % x = Data.StationID(northwest & strcmp(Data.Cruise, cruises.ID{7}));
        % [sum(strcmp(x, '11'))]

        % Cruise JR16003 contains 3 extra samples forming an eastward transect from
        % the identified station 11 -- label these as 11.1-11.3.
        cruiseID = 'JR16003';
        ind = strcmp(Data.Cruise, cruiseID) & position;
        ind = ind & cellfun(@(z) isempty(z), Data.StationID); % index all samples to label
        ind_ = false(ndat, 3);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,3) = ind & Data.Longitude == min(Data.Longitude(ind));
        ind_(:,2) = ind & ~(ind_(:,1) | ind_(:,3));
        for m = 1:3, Data.StationID(ind_(:,m)) = {[stationNames{:}, '.' num2str(m)]}; end


        % Cruise JR20120327 has 2 extra samples not associated to other stations --
        % label these as stations 12 and 13
        nstations = 2;
        whichCruise = strcmp(Data.Cruise, 'JR20120327');
        position = Data.Longitude < -41 & Data.Latitude > -54;
        ind = northwest & whichCruise & position; % indexes sample station from JR20141115
        % Find coordinates of these stations
        ind_ = false(ndat, nstations);
        ind_(:,1) = ind & Data.Longitude == max(Data.Longitude(ind));
        ind_(:,2) = ind & Data.Longitude == min(Data.Longitude(ind));
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % % Distance (km) between these stations
        % Dist = m_lldist(coords_(1,:), coords_(2,:));
        % radius = 0.5 * Dist; % choose a radius around sample stations -- samples within this radius matched to station

        % Define these as stations 12 and 13
        stationNames = cellstr(string(12:13));
        for m = 1:nstations, Data.StationID(ind_(:,m)) = stationNames(m); end

        % For all other cruises, identify whether this station was sampled then
        % set the station IDs to be consistent between cruises.
        % Find distances between all samples and the selected station
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Update the MetaData
        MetaData.StationID = cell(height(MetaData), 1);
        MetaData = movevars(MetaData, 'StationID', 'after', 'Station');
        for i = 1:height(MetaData)
            if ~strcmp(MetaData.Location{i}, 'northwest'), continue; end
            ind = strcmp(Data.Cruise, MetaData.Cruise{i}) & strcmp(Data.Label, MetaData.Label{i});
            MetaData.StationID(i) = unique(Data.StationID(ind));
        end

        % Replot the maps for the northwest cruises
        if plotMaps
            close(plotHandles.map2); close(plotHandles.map3); close(plotHandles.map4); close(plotHandles.map6); close(plotHandles.map7);
            for i = 1:height(cruises)
                if ~strcmp(cruises.Location{i}, 'northwest'), continue; end
                cruiseID = cruises.ID{i};
                md = MetaData(strcmp(MetaData.Cruise, cruiseID),:);
                lons = md.Longitude;
                lats = md.Latitude;
                labels = md.StationID;
                % Create South Georgia base map
                mapName = ['map' num2str(i)];
                figure
                plotHandles.(mapName) = gcf;
                %         assignin('base', mapName, figure)
                plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                    'XaxisLocation', XaxisLocation);
                % Extract map bounding coordinates
                [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0]);
                mapAx = gca;
                % Plot sample positions
                hold(mapAx, 'on')
                m_scatter(lons, lats, 18, [1 0 0], 'filled')
                % remove duplicate labels -- multiple samples at single station
                lt = table(lons, lats, labels);
                nlt = cell2mat(cellfun(@(z) str2num(z), lt.labels, 'UniformOutput', false));
                [~, o] = sort(nlt);
                lt = lt(o,:); nlt = nlt(o,:);
                ult = unique(nlt);
                for j = 1:length(ult)
                    k = nlt == ult(j);
                    x = lt(find(k,1),:);
                    lt(k,:) = []; nlt(k) = [];
                    lt = [lt; x]; nlt = [nlt; ult(j)];
                end
                m_text(lt.lons, lt.lats, lt.labels, 'FontSize', 8, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                titleLon = mean(lonMap);
                titleLat = 0.1 * diff(latMap) + max(latMap);
                switch ShowMapTitle, case true
                    m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' cruises.Season{i} ')'], ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                end
                hold(mapAx, 'off')
            end
        end

        % Station labels are now consistent between cruises sampling the northwest
        % region... now organise the southeast region

        % Use the 9 samples from cruise JR20150309 to set the station IDs.
        nstations = 9;
        stationNames = cellstr(string(1:nstations));

        cruiseID = 'JR20150309';
        whichCruise = strcmp(Data.Cruise, cruiseID);
        d = Data(whichCruise,:);
        samples = unique(d.Label);
        for i = 1:length(samples)
            sampleID = samples{i};
            whichSample = strcmp(d.Label, sampleID);
            d.StationID(whichSample) = repmat(stationNames(i), [sum(whichSample), 1]);
        end
        Data(whichCruise,:) = d;

        % position = Data.Longitude < -41 & Data.Latitude > -54;
        ind = southeast & whichCruise;
        % Find coordinates of these stations -- listed south -> north
        ind_ = false(ndat, nstations);
        for i = 1:nstations, ind_(:,i) = ind & strcmp(Data.StationID, num2str(i)); end
        ind_n = ind_ + 0; ind_n(ind_n == 0) = nan; % convert index from [0,1] logical -> [NaN,1] numeric
        lons = Data.Longitude .* ind_n;
        lats = Data.Latitude .* ind_n;
        coords_ = [mean(lons, 'omitnan'); mean(lats, 'omitnan')]; % 1 column per identified station

        % Get station IDs for other southeast cruises
        % Find distances between all samples and the selected stations
        lons = Data.Longitude; lats = Data.Latitude;
        lons_ = repmat(coords_(1,:), [2 * ndat + 1, 1]);
        lats_ = repmat(coords_(2,:), [2 * ndat + 1, 1]);
        lons_(2:2:end-1,:) = repmat(lons, [1, nstations]);
        lats_(2:2:end-1,:) = repmat(lats, [1, nstations]);
        allDist = nan(2 * ndat, nstations);
        for m = 1:nstations, allDist(:,m) = m_lldist(lons_(:,m), lats_(:,m)); end
        allDist = allDist(2:2:end,:); % distances between stations 1 & 2 and all other samples

        for m = 1:height(cruises)
            cruiseID = cruises.ID{m};
            ind = strcmp(Data.Cruise, cruiseID);
            d = allDist(ind,:); % distance between all cruise samples and the 2 identified stations
            nearStation = d < radius;
            if ~any(nearStation, 'all')
                continue
            else
                samples = unique(Data.Label(ind), 'stable');
                for s = 1:length(samples)
                    ind_ = ind & strcmp(Data.Label, samples{s});
                    d = allDist(ind_,:);
                    nearStation = d < radius;
                    if ~any(nearStation, 'all')
                        continue
                    else
                        whichStation = any(d == min(d, [], 2));
                        Data.StationID(ind_) = stationNames(whichStation);
                    end
                end
            end
        end

        % Update the MetaData
        for i = 1:height(MetaData)
            if ~strcmp(MetaData.Location{i}, 'southeast'), continue; end
            ind = strcmp(Data.Cruise, MetaData.Cruise{i}) & strcmp(Data.Label, MetaData.Label{i});
            MetaData.StationID(i) = unique(Data.StationID(ind));
        end

        % Replot the maps for the southeast cruises
        if plotMaps
            %     delete(map1); delete(map5)
            %     close('map1', 'map5')
            %     close(map1); close(map5)
            close(plotHandles.map1); close(plotHandles.map5);
            for i = 1:height(cruises)
                if ~strcmp(cruises.Location{i}, 'southeast'), continue; end
                cruiseID = cruises.ID{i};
                md = MetaData(strcmp(MetaData.Cruise, cruiseID),:);
                lons = md.Longitude;
                lats = md.Latitude;
                labels = md.StationID;
                % Create South Georgia base map
                mapName = ['map' num2str(i)];
                figure
                plotHandles.(mapName) = gcf;
                %         assignin('base', mapName, figure)
                plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                    'XaxisLocation', XaxisLocation);
                % Extract map bounding coordinates
                [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                    'redrawCoastline', false, 'edgecolour', [0 0 0]);
                mapAx = gca;
                % Plot sample positions
                hold(mapAx, 'on')
                m_scatter(lons, lats, 18, [1 0 0], 'filled')
                % remove duplicate labels -- multiple samples at single station
                lt = table(lons, lats, labels);
                nlt = cell2mat(cellfun(@(z) str2double(z), lt.labels, 'UniformOutput', false));
                [~, o] = sort(nlt);
                lt = lt(o,:); nlt = nlt(o,:);
                ult = unique(nlt);
                for j = 1:length(ult)
                    k = nlt == ult(j);
                    x = lt(find(k,1),:);
                    lt(k,:) = []; nlt(k) = [];
                    lt = [lt; x]; nlt = [nlt; ult(j)];
                end
                m_text(lt.lons, lt.lats, lt.labels, 'FontSize', 8, ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                titleLon = mean(lonMap);
                titleLat = 0.1 * diff(latMap) + max(latMap);
                switch ShowMapTitle, case true
                    m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' cruises.Season{i} ')'], ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                end
                hold(mapAx, 'off')
            end
        end


        %% Are samples depth profiles?

        % Not all of the samples are depth profiles -- many appear to be
        % measurements at/around a single depth, often near the seafloor.
        % Determine which sample are depth profiles and index them.

        % This can be automated using coefficient of variation of measurement depth
        % -- very small CoV => measurements cluster around single depth.

        % CoV_lim = 0.05; % CoV >= CoV_lim => depth profile, CoV < CoV_lim => sample around single depth. CoV_lim is arbitrary choice...
        Data.DepthCoV = nan(height(Data),1);
        Data.DepthProfile = false(height(Data),1);
        MetaData.DepthCoV = false(height(MetaData), 1);
        MetaData.DepthProfile = false(height(MetaData), 1);

        for i = 1:height(cruises)
            whichCruise = cruises.ID{i};
            indi = strcmp(Data.Cruise, whichCruise); % index cruise
            allStations = unique(Data.StationID(indi), 'stable');
            for j = 1:length(allStations)
                whichStation = allStations{j};
                indj = indi & strcmp(Data.StationID, whichStation);  % index cruise & station
                allSamples = unique(Data.Label(indj));
                for k = 1:length(allSamples)
                    whichSample = allSamples{k};
                    indk = indj & strcmp(Data.Label, whichSample);
                    indm = strcmp(MetaData.Cruise, whichCruise) & strcmp(MetaData.StationID, whichStation) & strcmp(MetaData.Label, whichSample);
                    CoV = std(Data.Depth(indk), 'omitnan') / mean(Data.Depth(indk), 'omitnan');
                    Data.DepthCoV(indk) = CoV;
                    MetaData.DepthCoV(indm) = CoV;
                    if CoV >= CoV_lim
                        Data.DepthProfile(indk) = true;
                        MetaData.DepthProfile(indm) = true;
                    end
                    %             if i ~= 1 || j ~= 1 || k ~= 1
                    %                 profileSummary.Cruise = [profileSummary.Cruise; whichCruise];
                    %                 profileSummary.StationID = [profileSummary.StationID; whichStation];
                    %                 profileSummary.Label = [profileSummary.Label; whichSample];
                    %                 profileSummary.DepthCoV = [profileSummary.DepthCoV; CoV];
                    %                 profileSummary.DepthProfile = [profileSummary.DepthProfile; CoV >= CoV_lim];
                    %             else
                    %                 profileSummary.Cruise = {whichCruise};
                    %                 profileSummary.StationID = {whichStation};
                    %                 profileSummary.Label = {whichSample};
                    %                 profileSummary.DepthCoV = CoV;
                    %                 profileSummary.DepthProfile = CoV >= CoV_lim;
                    %             end
                end
            end
        end

        % profileSummary = struct2table(profileSummary);

        t = sprintf('%u%%', round(100 * sum(MetaData.DepthProfile) / height(MetaData)));
        t = sprintf('Of all CTD samples, %s', t);
        fprintf('\n\n'); fprintf(1, '%s are depth profiles.\nThe remaining samples are measurements around a single depth.\n(Automatic determination of depth profiles depends on CoV_lim,\nwhich is a somewhat arbitrary choice so maybe try ajusting\nthis parameter... lower values include more samples.)', t); fprintf('\n\n')


        %% Flag spurious points

        % Although the quality control flags suggest all measurements are valid,
        % plotting demonstrates lots of spurious measurements that mostly appear
        % in the first few data points. This suggests that sensor readings are
        % unreliable when CTD is first deployed, so these measurements should be
        % omitted.
        % Focus on depth profile samples -- ignoring all others for now.

        % Each sample is a series of successive measurements seprated by some,
        % hopefully uniform, time increment.
        % Plot measurement depth against measurement number -- separate multipanel
        % plot for each cruise.

        if plotSpuriousPoints
            for i = 1:height(cruises)
                whichCruise = cruises.ID{i};
                ps = MetaData(strcmp(MetaData.Cruise, whichCruise),:);
                if ~any(ps.DepthProfile), continue; end
                np = sum(ps.DepthProfile);
                nrows = ceil(np ^ 0.5); ncols = ceil(np / nrows);
                testPlotHandles.(['plot' num2str(i)]) = figure;
                set(gcf, {'Unit','Position'}, {'inches', [0 0 4*ncols 3*nrows]})
                j = 0;
                for k = 1:height(ps)
                    if ~ps.DepthProfile(k), continue; end
                    j = j + 1;
                    subplot(nrows, ncols, j)
                    whichStation = ps.StationID{k};
                    whichSample = ps.Label{k};
                    ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation) & strcmp(Data.Label, whichSample);
                    Dat = Data(ind,:);
                    sampleNumber = 1:height(Dat);
                    scatter(sampleNumber, -Dat.Depth)
                    xlabel('Measurement'); ylabel('Depth (m)'); title([whichCruise ' - station ' whichStation ' - ' whichSample], 'FontWeight', 'normal')
                end
            end
        end

        % Clearly, many depth profiles contain spurious values, mostly appearing as
        % the initial measurements.
        % Points forming lines with negative gradient correspond to good
        % measurements along a depth profile.

        % The following is an ad-hoc method to catch spurious points... it's imperfect

        Data.SpuriousMeasure = false(height(Data),1);

        % This is a decent method, but it doesn't quite catch all the spurious
        % points... Just manually find the spurious points, using the following
        % plots as a guide...

        for i = 1:height(cruises)
            whichCruise = cruises.ID{i};
            indi = strcmp(Data.Cruise, whichCruise); % index cruise
            allStations = cellstr(string(unique(cellfun(@(z) str2double(z), Data.StationID(indi)))));
            for j = 1:length(allStations)
                whichStation = allStations{j};
                indj = indi & strcmp(Data.StationID, whichStation); % index cruise & station
                allSamples = unique(Data.Label(indj));
                for k = 1:length(allSamples)
                    whichSample = allSamples{k};
                    indk = indj & strcmp(Data.Label, whichSample); % index cruise & station & sample
                    isDepthProfile = MetaData.DepthProfile(strcmp(MetaData.Cruise, whichCruise) & strcmp(MetaData.StationID, whichStation) & strcmp(MetaData.Label, whichSample));
                    if ~isDepthProfile, continue; end
                    Dat = Data(indk,:);
                    np = height(Dat); % number of measurements
                    mp = 1:np;
                    badPoints = false(np,1); % index spurious data points
                    %             figure
                    %             scatter(1:sum(indk), -Dat.Depth); xlabel('Measurement'); ylabel('Depth (m)')

                    % Are the final measurements OK? The last point should be the
                    % deepest. If it's not then there's a problem...
                    deepestPoint = Dat.Depth == max(Dat.Depth);
                    if ~deepestPoint(end)
                        badPoints = badPoints | logical(cumsum([0; deepestPoint(1:end-1)]));
                    end
                    nbp = sum(badPoints);

                    %             hold on
                    %             scatter(mp(badPoints), -Dat.Depth(badPoints), 'MarkerEdgeColor', [1 0 0])
                    %             hold off
                    % Spurious final measurements have (hopefully) been identified,
                    % but note that this method will fail if the deepest samples
                    % appear as spurious initial measurements... it looks OK for
                    % the BODC South Georgia.
                    % Now remove spurious initial measurements.
                    d = flip(Dat.Depth);
                    badPoints = flip(badPoints);
                    %             scatter(mp, d) % isolate points on the linear slope of negative gradient
                    %             hold on
                    %             scatter(mp(badPoints), d(badPoints), 'MarkerEdgeColor', [1 0 0])
                    %             hold off
                    d_ = d(~badPoints); % omit any spurious final points
                    dd = diff(d_); dd = [dd(1); dd];
                    ma = ceil(movstdRange * np);
                    mm = movmean(dd, [ma, 0]);
                    ms = movstd(dd, [ma, 0]);
                    mm = [mm(1); mm(1:end-1)]; % shift by one to get expected values
                    ms = [ms(1); ms(1:end-1)];
                    %             ms(1:nbp+1) = ms(nbp+2); % remove zeros

                    cutOffMetric = zeros(np, leeway+1);
                    cutOffMetric(nbp+1:end,1) = abs(ms ./ dd); % .^ 2;
                    for l = 1:leeway
                        cutOffMetric(nbp+1:end,l+1) = abs([ms(1:end-l) ./ dd(l+1:end); zeros(l, 1)]);
                    end
                    cutOffMetric(1:ma,:) = 0; % avoid false positives at final measurements
                    badPoints = badPoints | cumsum(all(cutOffMetric > cutOffValue, 2)) > 0;
                    badPoints = flip(badPoints);
                    %             scatter(mp, -Dat.Depth); xlabel('Measurement'); ylabel('Depth (m)')
                    %             hold on
                    %             scatter(mp(badPoints), -Dat.Depth(badPoints), 'MarkerEdgeColor', [1 0 0])
                    %             hold off
                    Data.SpuriousMeasure(indk) = badPoints;
                end
            end
        end

        if plotSpuriousPoints
            % Replot to view spurious measurements
            fields = fieldnames(testPlotHandles);
            for i = 1:length(fields), close(testPlotHandles.(fields{i})); end
            clearvars testPlotHandles

            for i = 1:height(cruises)
                whichCruise = cruises.ID{i};
                ps = MetaData(strcmp(MetaData.Cruise, whichCruise),:);
                if ~any(ps.DepthProfile), continue; end
                np = sum(ps.DepthProfile);
                nrows = ceil(np ^ 0.5); ncols = ceil(np / nrows);
                testPlotHandles.(['plot' num2str(i)]) = figure;
                set(gcf, {'Unit','Position'}, {'inches', [0 0 4*ncols 3*nrows]})
                j = 0;
                for k = 1:height(ps)
                    if ~ps.DepthProfile(k), continue; end
                    j = j + 1;
                    subplot(nrows, ncols, j)
                    whichStation = ps.StationID{k};
                    whichSample = ps.Label{k};
                    ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation) & strcmp(Data.Label, whichSample);
                    Dat = Data(ind,:);
                    sampleNumber = 1:height(Dat);
                    scatter(sampleNumber, -Dat.Depth)
                    xlabel('Measurement'); ylabel('Depth (m)'); title([whichCruise ' - station ' whichStation ' - ' whichSample], 'FontWeight', 'normal')
                    hold on
                    spuriousPoints = Dat.SpuriousMeasure;
                    scatter(sampleNumber(spuriousPoints), -Dat.Depth(spuriousPoints), 'MarkerEdgeColor', [1 0 0])
                    hold off
                end
            end
        end

        % Now, with help from these plots, manually specify the spurious points...
        switch source
            case 'BODC'
                whichCruise = 'JR20120327'; % no spurious points
                Data.SpuriousMeasure(strcmp(Data.Cruise, whichCruise)) = false;
                whichCruise = 'JR20131112'; % the above filtering method worked OK
                whichCruise = 'JR20141115'; % stations 3 & 4 need fixed
                whichStation = '3';
                ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation);
                x = false(sum(ind), 1); x(1:16) = true;
                Data.SpuriousMeasure(ind) = x;
                whichStation = '4';
                ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation);
                x = false(sum(ind), 1); x(1:57) = true;
                Data.SpuriousMeasure(ind) = x;
                whichCruise = 'JR15002'; % the above filtering method worked OK
                whichCruise = 'JR16003'; % the above filtering method worked OK
        end


        if plotSpuriousPoints

            % Replot to view spurious measurements
            fields = fieldnames(testPlotHandles);
            for i = 1:length(fields), close(testPlotHandles.(fields{i})); end
            clearvars testPlotHandles

            for i = 1:height(cruises)
                whichCruise = cruises.ID{i};
                ps = MetaData(strcmp(MetaData.Cruise, whichCruise),:);
                if ~any(ps.DepthProfile), continue; end
                np = sum(ps.DepthProfile);
                nrows = ceil(np ^ 0.5); ncols = ceil(np / nrows);
                testPlotHandles.(['plot' num2str(i)]) = figure;
                set(gcf, {'Unit','Position'}, {'inches', [0 0 4*ncols 3*nrows]})
                j = 0;
                for k = 1:height(ps)
                    if ~ps.DepthProfile(k), continue; end
                    j = j + 1;
                    subplot(nrows, ncols, j)
                    whichStation = ps.StationID{k};
                    whichSample = ps.Label{k};
                    ind = strcmp(Data.Cruise, whichCruise) & strcmp(Data.StationID, whichStation) & strcmp(Data.Label, whichSample);
                    Dat = Data(ind,:);
                    sampleNumber = 1:height(Dat);
                    scatter(sampleNumber, -Dat.Depth)
                    xlabel('Measurement'); ylabel('Depth (m)'); title([whichCruise ' - station ' whichStation ' - ' whichSample], 'FontWeight', 'normal')
                    hold on
                    spuriousPoints = Dat.SpuriousMeasure;
                    scatter(sampleNumber(spuriousPoints), -Dat.Depth(spuriousPoints), 'MarkerEdgeColor', [1 0 0])
                    hold off
                end
                sgtitle('Depth profiles with spurious points highlighted')
            end

        end


        %% Tidy up

        % Check for and omit any negative (above surface) depths
        Data = Data(Data.Depth > 0,:);

        % Remove CTD casts that are not depth profiles?
        switch onlyDepthProfiles, case true
            Data = Data(Data.DepthProfile,:);
            MetaData = MetaData(MetaData.DepthProfile,:);
            fprintf('\n\n'); fprintf('The data have been filtered to only include CTD depth profiles.\nAll other CTD casts are excluded and not shown on maps.'); fprintf('\n\n')

            % Replot maps to show only depth profiles
            if plotMaps
                fields = fieldnames(plotHandles);
                for i = 1:length(fields), close(plotHandles.(fields{i})); end
                clearvars plotHandles

                % Create a separate map for each cruise
                ucruises = unique(MetaData.Cruise, 'stable');
                ncruises = length(ucruises);
                for i = 1:ncruises
                    cruiseID = ucruises{i};
                    seasonID = cruises.Season{strcmp(cruises.ID, cruiseID)};
                    Dat = Data(strcmp(Data.Cruise, cruiseID),:);
                    mDat = MetaData(strcmp(MetaData.Cruise, cruiseID),:);
                    labels = mDat.StationID;
                    %             samples = mDat.Label;
                    nsamples = length(labels);
                    lons = mDat.Longitude; lats = mDat.Latitude;
                    %             labels = cellfun(@(z) num2str(str2double(z(4:end))), mDat.Label, 'UniformOutput', false);

                    % Create South Georgia base map
                    mapName = ['map' num2str(i)];
                    figure
                    plotHandles.(mapName) = gcf;
                    %             assignin('base', mapName, figure)
                    plotBaseMap('South Georgia', 'createMap', true, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0], ...
                        'XaxisLocation', XaxisLocation);
                    % Extract map bounding coordinates
                    [lonMap, latMap] = plotBaseMap('South Georgia', 'createMap', false, 'coordsTable', coordsTable, ...
                        'redrawCoastline', false, 'edgecolour', [0 0 0]);
                    mapAx = gca;
                    % Plot sample positions
                    hold(mapAx, 'on')
                    m_scatter(lons, lats, 18, [1 0 0], 'filled')
                    m_text(lons, lats, labels, 'FontSize', 8, ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                    titleLon = mean(lonMap);
                    titleLat = 0.1 * diff(latMap) + max(latMap);
                    switch ShowMapTitle, case true
                        m_text(titleLon, titleLat, ['CTD casts: cruise ' cruiseID, ' (' seasonID ')'], ...
                            'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'FontSize', 12)
                    end
                    hold(mapAx, 'off')
                end
                clearvars counter Dat samples nsamples lons lats labels lonMap latMap mapAx titleLon titleLat i j k l *ID
            end
        end


        % Remove spurious points?
        switch removeSpuriousPoints
            case true
                Data = Data(~Data.SpuriousMeasure,:);
            case false
                disp('Spurious measurements may remain in the data. Examine "Data.SpuriousMeasure" to check.')
        end




end


