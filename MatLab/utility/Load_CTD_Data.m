function [Data, MetaData, cruises] = Load_CTD_Data(project, source, region, varargin)
% Load and organise CTD data sampled from Region and downloaded from Source

extractVarargin(varargin)

if ~exist('printCruiseStructure', 'var'), printCruiseStructure = false; end % true => structure containing cruise metadata is printed on screen
if ~exist('printFileStructure', 'var'), printFileStructure = false; end % true => structure containing data file names is displayed on screen
if ~exist('displayData', 'var'), displayData = true; end

% Set top level directory
thisFile = which('Load_CTD_Data');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);

% File directory structure depends on data source so use switch-case to
% compartmentalise code.
switch source
    case 'BODC'
        %% Set directories and basic structs
        % Data are stored in separate directories for each season -- find available
        % seasons
        df = dir(fullfile(baseDirectory, ...
            'data', 'CTD', source, region));
        df = df(~ismember({df.name},{'.','..'}));

        seasons.nseasons = length({df.name});
        seasons.season = {df.name};
        seasons.fieldName = strcat('season', strrep(seasons.season, '-', '_')); % labels for structs
        clearvars df; fprintf(1, '\n\nData are available from seasons: %s', [seasons.season{1}, sprintf(', %s', seasons.season{2:end})]); fprintf(1, '\n\n\n')

        % Find cruise IDs from directory names -- data from each cruise is stored
        % in separate directory
        for i = 1:seasons.nseasons
            if i == 1, cruises.ncruises = 0; end
            df = fullfile(baseDirectory, ...
                'data', 'CTD', source, region, seasons.season{i});
            d = dir(df);
            d = d(~ismember({d.name},{'.','..'}));
            d_ = struct2cell(d);
            cruises.ncruises = cruises.ncruises + size(d_, 2);
            cruises.(seasons.fieldName{i}).ncruises = size(d_, 2);
            cruises.(seasons.fieldName{i}).name = d_(strcmp(fieldnames(d), 'name'),:);
            cruises.(seasons.fieldName{i}).ID = strrep(cruises.(seasons.fieldName{i}).name, 'cruise ','');
            clearvars df d d_ i
        end
        switch printCruiseStructure, case true
            fprintf('\nAll seasons...\n\n')
            disp(cruises)
            fprintf('First listed season...\n\n')
            disp(cruises.(seasons.fieldName{1})); clearvars printCruiseStructure;
            otherwise, clearvars printCruiseStructure
        end

        % Specify data directories
        for i = 1:seasons.nseasons
            df = fullfile(baseDirectory, ...
                'data', 'CTD', source, region, seasons.season{i});
            ncruises = cruises.(seasons.fieldName{i}).ncruises;
            cruiseNames = cruises.(seasons.fieldName{i}).name;
            cruiseIDs = cruises.(seasons.fieldName{i}).ID;
            for j = 1:ncruises
                dfj = fullfile(df, cruiseNames{j});
                dirData.(seasons.fieldName{i}).(cruiseIDs{j}) = dfj;
                addpath(dfj)
            end
            clearvars df dfj ncruises cruiseNames cruiseIDs i j
        end


        %% Load CTD data

        % Some data files appear to be samples from a single depth rather than
        % profiles... this may take a while to organise...

        % NetCDF formatting is not identical across seasons, so take care loading

        fileType = 'nc'; % data stored as NetCDF files
        % Find all data files of fileType stored in dirData
        for i = 1:seasons.nseasons
            seasonID = seasons.fieldName{i};
            dc = cruises.(seasonID);
            for j = 1:dc.ncruises
                cruiseID = dc.ID{j};
                df = dirData.(seasonID).(cruiseID);
                files = dir([df '/*.' fileType]); % list all (.fileType) files
                nfiles = length(files);
                fileNames.(seasonID).(cruiseID) = cell(nfiles, 1);
                for k = 1:nfiles, fileNames.(seasonID).(cruiseID){k} = files(k).name; end
            end
            clearvars seasonID dc cruiseID df files nfiles i j k
        end
        switch printFileStructure, case true
            fprintf('\nFile name structure...\n\n')
            disp(fileNames)
            fprintf('within a season...\n\n')
            disp(fileNames.(seasons.fieldName{1}))
            fprintf('within a cruise\n\n')
            disp(fileNames.(seasons.fieldName{1}).(cruises.(seasons.fieldName{1}).ID{1})); clearvars printFileStructure;
            otherwise, clearvars printFileStructure
        end

        % Extract data fields from saved NetCDF files
        printDataStructure = displayData; % true => data structure outline is shown on screen after compilation
        for i = 1:seasons.nseasons
            if i == 1, fprintf('\n\nLoading NetCDF data...\n\n'); end
            seasonID = seasons.fieldName{i};
            dc = cruises.(seasonID);
            for j = 1:dc.ncruises
                cruiseID = dc.ID{j};
                fn = fileNames.(seasonID).(cruiseID);
                if isempty(fn), continue; end % if no files in directory then skip
                for k = 1:length(fn)
                    fileName = fn{k};
                    datInfo = ncinfo(fileName);
                    variableNames = {datInfo.Variables.Name}; % NOTE: variable names are not consistent across seasons
                    %         ncdisp(fileName)
                    Dat.Cruise = cruiseID;
                    Dat.Station = ncread(fileName, 'SDN_STATION');
                    Dat.Station = cellstr(Dat.Station');
                    Dat.Station = Dat.Station{1};
                    Dat.Label = ['CTD' num2str(k)];
                    Dat.Time = ncread(fileName, 'TIME');
                    Dat.Longitude = ncread(fileName, 'LONGITUDE');
                    Dat.Latitude = ncread(fileName, 'LATITUDE');
                    Dat.SeaFloorDepth = ncread(fileName, 'SDN_BOT_DEPTH'); % unit: m. Sea-floor depth (below instantaneous sea level) {bathymetric depth} in the water body
                    Dat.HeightAboveSeaFloor = ncread(fileName, 'AHSFZZ01'); % unit: m. Height (spatial coordinate) relative to bed surface in the water body
                    Dat.HeightAboveSeaFloorFlag = ncread(fileName, 'AHSFZZ01_SEADATANET_QC');
                    Dat.Depth = Dat.SeaFloorDepth - Dat.HeightAboveSeaFloor;
                    Dat.DepthFlag = Dat.HeightAboveSeaFloorFlag;
                    Dat.PressureSeawater = ncread(fileName, 'PRES'); % unit: dbar. Pressure (spatial coordinate) exerted by the water body by profiling pressure sensor and correction to read zero at sea level
                    Dat.PressureSeawaterFlag = ncread(fileName, 'PRES_SEADATANET_QC');
                    possibleNames = {'CNDCST01', 'CNCLCCI1'}; % conductivity labels differ between seasons
                    varName = possibleNames{contains(possibleNames, variableNames)};
                    Dat.Conductivity = ncread(fileName, varName); % unit: S/m
                    Dat.ConductivityFlag = ncread(fileName, [varName '_SEADATANET_QC']);
                    possibleNames = {'CPHLPR01', 'CPHLPM01'}; % chlorophyll labels differ between seasons
                    varName = possibleNames{contains(possibleNames, variableNames)};
                    Dat.Chlorophyll = ncread(fileName, varName); % unit: mg chl-a/m^3
                    Dat.ChlorophyllFlag = ncread(fileName, [varName '_SEADATANET_QC']);
                    Dat.TemperaturePotential = ncread(fileName, 'POTMCV01'); % unit: degree C. Potential temperature of the water body by computation using UNESCO 1983 algorithm
                    Dat.TemperaturePotentialFlag = ncread(fileName, 'POTMCV01_SEADATANET_QC');
                    possibleNames = {'PSALST01', 'PSALCU01'}; % salinity labels differ between seasons
                    varName = possibleNames{contains(possibleNames, variableNames)};
                    Dat.SalinityPractical = ncread(fileName, varName); % unit: dimensionless. Practical salinity of the water body by CTD and computation using UNESCO 1983 algorithm
                    Dat.SalinityPracticalFlag = ncread(fileName, [varName '_SEADATANET_QC']);
                    Dat.Density = ncread(fileName, 'SIGTPR01'); % unit: kg/m^3. Sigma-theta of the water body by CTD and computation from salinity and potential temperature using UNESCO algorithm
                    Dat.DensityFlag = ncread(fileName, 'SIGTPR01_SEADATANET_QC');
                    possibleNames = {'TEMPST01', 'TEMPCC01'}; % temperature labels differ between seasons
                    varName = possibleNames{contains(possibleNames, variableNames)};
                    Dat.Temperature = ncread(fileName, varName); % unit: degree C. Temperature of the water body by CTD or STD
                    Dat.TemperatureFlag = ncread(fileName, [varName '_SEADATANET_QC']);
                    Data.(seasonID).(cruiseID).(['CTD' num2str(k)]) = Dat; % store all data

                    nDat = fieldnames(Dat);
                    remove = structfun(@(z) size(z, 1) > 1, Dat);
                    Dat = rmfield(Dat, nDat(remove));
                    Dat.Season = seasons.season{i};
                    Dat = orderfields(Dat, [{'Season'}; nDat(~remove)]);
                    fields = fieldnames(Dat);
                    for m = 1:length(fields)
                        if ischar(Dat.(fields{m})), Dat.(fields{m}) = {Dat.(fields{m})}; end
                    end
                    Dat = struct2table(Dat);
                    if ~(i == 1 && j == 1 && k == 1)
                        MetaData = [MetaData; Dat]; % store (scalar) metadata separately
                    else
                        MetaData = Dat;
                    end
                    clearvars Dat nDat remove datInfo possibleNames variableNames varName nDat m fields
                end
            end
            %     fprintf(1, '%s / 100', num2str(round(100 * i / seasons.nseasons, 3, 'significant'))); fprintf('\n')
            clearvars dc fn fileName seasonID cruiseID i j k
        end
        switch printDataStructure, case true
            fprintf('\n\nData are stored as nested struct with 3 levels: season, cruise, and CTD cast.\n\n1: season\n\n');
            disp(Data); fprintf('\n2: cruise\n\n')
            disp(Data.(seasons.fieldName{1})); fprintf('\n3: CTD cast\n\n')
            disp(Data.(seasons.fieldName{1}).(cruises.(seasons.fieldName{1}).ID{1})); fprintf('\nMeasurements from each CTD cast stored separately within third level of nested struct.\n\n')
            disp(Data.(seasons.fieldName{1}).(cruises.(seasons.fieldName{1}).ID{1}).CTD1); clearvars printDataStructure; otherwise, clearvars printDataStructure
        end

    case 'BAS central storage'
        %% Set directories and basic structs
        if isempty(region)
            df = dir(fullfile(baseDirectory, ...
                'data', 'CTD', source, region));
        else
            df = dir(fullfile(baseDirectory, ...
                'data', 'CTD', source));
        end
        df = df(~ismember({df.name},{'.','..'}));

        % For now I'm only using data from one cruise so the code will
        % simple relative to the above BODC case.
        % The cruise is JR82 in 2003.
        cruise = 'JR82';
        season = 2003;
        dirCTDtimes = fullfile(df.folder, df.name);
        dirCTD = fullfile(df.folder, df.name, 'matlab');


        %% Load CTD data
        Data = load('jr82data_nospikes.mat');
        CTDtimes = readcell('sample times.txt', 'Delimiter', ' ');

%         % These data are in terms of pressure and exclude depth.
%         % Calculate a depth field using the hydrostatic water pressure
%         % formula: pressure = density * gravity acceleratoin * depth
%         rho = 1023.6; % saltwater density (at surface??)
%         g = 9.81;
%         P = Data.jr82press; % pressure in decibars -- I think!
%         P = P * 1e4; % convert decibars to Pascals
%         d = P ./ (rho * g); % this shoud be improved... the density is not independent of depth
%         Data.depth = d;


%         Data.jr82datenum = nan(size(Data.jr82stn));
        CTDevent = cell2mat(CTDtimes(2:end,1));
        
        for i = length(Data.jr82stn):-1:1
            stn = Data.jr82stn(i);
            j = [false; CTDevent == stn];
            time = CTDtimes{j,3};
            Data.jr82date(i,:) = time;
            Data.jr82datenum(i,:) = datenum(time);
        end

        % Remove the jr82 prefix from field name
        fields = fieldnames(Data);
        for i = 1:length(fields)
            if strcmp(fields{i}(1:4), 'jr82')
                Data = renameStructField(Data, fields{i}, fields{i}(5:end));
            end
        end

        MetaData = [];
        cruises = [];
end
