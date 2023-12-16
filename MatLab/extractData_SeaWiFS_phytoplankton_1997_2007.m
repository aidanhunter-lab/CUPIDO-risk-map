% Extract SeaWiFS global plankton data from netCDF files.
% Each file is about 0.25 GB and there are twelve -- one per month

% Data downloaded from:
% https://data.ceda.ac.uk/neodc/nceo-carbon/data/ST6-ocean_biogeochemistry/Global-PSC-Climatologies/monthly-10yr-climatology/

% If the netCDF files have already been downloaded then they should be
% saved under default file names in directory:
% ../CUPIDO-risk-map/data/chlorophyll/SeaWiFS/

% If the .nc files are already saved then this script will run quickly.
% Otherwise the files will be (temporarily) download, the data extracted 
% and saved as .mat, then the .nc files can be deleted -- this takes longer
% and requires an internet connection, but there's no need to store large
% data files locally.
addpath(genpath(fileparts(which('extractData_SeaWiFS_phytoplankton_1997_2007'))))
dirBase = fileparts(fileparts(which('extractData_SeaWiFS_phytoplankton_1997_2007')));
dirData = fullfile(dirBase, 'data', 'chlorophyll', 'SeaWiFS');
addpath(dirData)
% Unique segment of data download links
datalink_ = 'https://dap.ceda.ac.uk/neodc/nceo-carbon/data/ST6-ocean_biogeochemistry/Global-PSC-Climatologies/monthly-10yr-climatology//SeaWiFS_GLOBAL_19970901_Phytoplankton-Size-Class-10yr-Climatology-';
% Unique segment of default file names
filename_ = 'SeaWiFS_GLOBAL_19970901_Phytoplankton-Size-Class-10yr-Climatology-';
% Separate data set for each month
months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
% Store large .nc file locally?
saveNetCDF = true;

% https://dap.ceda.ac.uk/neodc/nceo-carbon/data/ST6-ocean_biogeochemistry/Global-PSC-Climatologies/monthly-10yr-climatology/

% Choose coordinate range to extract -- this covers whole Southern Ocean
lonrange = [-180, 180];
latrange = [-90, -45];

% Load, extract, and store data
% (This script could be greatly speeded up by using the OPenNDAP data
% filtering before downloading. See how the data URL changes when filtering
% parameters are selected in the online data access form. I did that in
% some other download scripts in this code repo, but must have created this
% script before learning how to do it! Maybe worth rejigging this script to
% improve download speeds?)

dataIsStored = false(1, length(months));
for m = 1:length(months)
    month = months{m};
    filename = [filename_, month, '.nc'];
    filepath = fullfile(dirData, filename);
    dataIsStored(m) = exist(filename, 'file') == 2;
    switch dataIsStored(m), case false
        datalink = [datalink_, month, '.nc'];
        disp(['Downloading from: ', datalink])
        websave(filepath, datalink); % download data
    end
    d = ncinfo(filename); % data details
    dims = {d.Dimensions.Name; d.Dimensions.Length}; % number of unique lat-longs
    nlon = dims{2,strcmp(dims(1,:), 'longitude')};
    nlat = dims{2,strcmp(dims(1,:), 'latitude')};
    dat.month = repmat(single(m), [nlon, nlat]);
    variables = {d.Variables.Name};
    for i = 1:length(variables)
        dat.(variables{i}) = ncread(filename, variables{i});
        if all(flip(size(dat.(variables{i}))) == size(dat.month)), dat.(variable{1}) = permute(dat.(variables{i}), [2, 1]); end
    end
    ind = lonrange(1) <= dat.longitude & dat.longitude <= lonrange(2) & ...
        latrange(1) <= dat.latitude & dat.latitude <= latrange(2);
    x = any(ind, 2);
    y = any(ind, 1);
    if length(x) == nlon, nlon_ = sum(x); else, nlon_ = sum(y); end
    if length(x) == nlat, nlat_ = sum(x); else, nlat_ = sum(y); end
    dat = structfun(@(z) reshape(z(ind), [nlon_, nlat_]), dat, 'UniformOutput', false);
    if m ~= 1
        fields = fieldnames(dat);
        for i = 1:length(fields)
            Dat.(fields{i}) = cat(3, Dat.(fields{i}), dat.(fields{i}));
        end
    else
        Dat = dat;
    end
    switch saveNetCDF, case false
        delete(filepath)
    end
    clear dat
end

wasDataDownloaded = any(~dataIsStored);
switch wasDataDownloaded, case true
    downloadTime = {date};
    note = cell2table(downloadTime);
    writetable(note, fullfile(dirData, 'time of data download.txt'))
end

disp(Dat)

% Save the data as a MatLab structure (.mat) and/or as a table (.csv or
% .txt).
% File sizes are large so only save data in native MatLab (.mat) format
ext = '.mat';
outputFile = ['SeaWiFS_Phytoplankton-Size-Class-1997-2007-Southern-Ocean', ext];
outputPath = fullfile(dirData, outputFile);

switch ext
    case '.mat'
        save(outputPath, 'Dat')
    case {'.txt', '.csv'}
        Dat = structfun(@(z) z(:), Dat, 'UniformOutput', false); % table variables change in order: longitude, latitude, month
        Dat = struct2table(Dat);
        writetable(Dat, outputPath)
end


%% Subset regions from Southern Ocean & save separately as NetCDF or csv
% (This code, below, is not necessary, it just allows saving relatively
% small NetCDF files for some specific regions listed in 
% 'regional bounding coordinates.csv'.

coordsTable = readtable('regional bounding coordinates.csv');

Dat_ = Dat;
nlocations = size(coordsTable, 1);
locations = cell(nlocations, 1);
for i = 1:size(coordsTable, 1)
    location = coordsTable.Location{i};
    location = strrep(location, ',', '');
    location = strrep(location, ' ', '_');
    locations{i} = location;
    lon = [coordsTable.Longitude_min(i), coordsTable.Longitude_max(i)];
    lat = [coordsTable.Latitude_min(i), coordsTable.Latitude_max(i)];    
    ind = lat(1) <= Dat.latitude & Dat.latitude <= lat(2) & ...
        lon(1) <= Dat.longitude & Dat.longitude <= lon(2);    
    nlon = max(max(sum(ind)));
    nlat = max(max(sum(ind,2)));
    Dat.(location) = structfun(@(z) reshape(z(ind), [nlon, nlat, 12]), Dat_, 'UniformOutput', false);
end
clearvars Dat_

% Save regional data as NetCDF files
outputFile_ = 'SeaWiFS_Phytoplankton-Size-Class-1997-2007-';
ext = '.nc';

for i = 1:size(coordsTable, 1)
    location = locations{i};
    outputFile = [outputFile_, location, ext];
    outputPath = fullfile(dirData, outputFile);
    dat = Dat.(location);
    
    lon = squeeze(dat.longitude(:,1,1));
    lat = squeeze(dat.latitude(1,:,1));
    mon = squeeze(dat.month(1,1,:));
    
    nccreate(outputPath, 'longitude', 'Dimensions', {'longitude', length(lon), 'latitude', length(lat), 'month', length(mon)}, 'Datatype', 'single')
    nccreate(outputPath, 'latitude', 'Dimensions', {'longitude', 'latitude', 'month'}, 'Datatype', 'single')
    nccreate(outputPath, 'month', 'Dimensions', {'longitude', 'latitude', 'month'}, 'Datatype', 'single')
    nccreate(outputPath, 'Tchl', 'Dimensions', {'longitude', 'latitude', 'month'}, 'Datatype', 'single')
    nccreate(outputPath, 'Pico_percent_Tchl', 'Dimensions', {'longitude', 'latitude', 'month'}, 'Datatype', 'single')
    nccreate(outputPath, 'Nano_percent_Tchl', 'Dimensions', {'longitude', 'latitude', 'month'}, 'Datatype', 'single')
    nccreate(outputPath, 'Micro_percent_Tchl', 'Dimensions', {'longitude', 'latitude', 'month'}, 'Datatype', 'single')
    
    schema.all.Name = '/';
    attributeNames = {'Title', 'Institution', 'Source', 'History', 'References', 'Data subset'};
    attributeValues = {'Global data on Phytoplankton Size Structure as observed from the SeaWiFS satellite',...
        'Plymouth Marine Laboratory', ...
        'NASA Ocean Color website http://oceancolor.gsfc.nasa.gov/', ...
        'April 2012: Level 3 mapped global 9km Chlorophyll-a data downloaded from NASA Ocean Color website and the model of Brewin et al. (2010) applied on a pixel-by-pixel basis to retrieve size structure', ...
        'Brewin, R.J.W., Sathyendranath, S., Hirata, T., Lavender, S., Barciela, R.M. & Hardman- Mountford, N.J. (2010). A three-component model of phytoplankton size class for the Atlantic Ocean. Ecological Modelling, 221, 1472−1483.', ...
        location};
    schema.all.Attributes = struct('Name', attributeNames, 'Value', attributeValues);

    ncwriteschema(outputPath, schema.all)
    
    ncwrite(outputPath, 'longitude', dat.longitude)
    ncwrite(outputPath, 'latitude', dat.latitude)
    ncwrite(outputPath, 'month', dat.month)
    ncwrite(outputPath, 'Tchl', dat.Tchl)
    ncwrite(outputPath, 'Pico_percent_Tchl', dat.Pico_percent_Tchl)
    ncwrite(outputPath, 'Nano_percent_Tchl', dat.Nano_percent_Tchl)
    ncwrite(outputPath, 'Micro_percent_Tchl', dat.Micro_percent_Tchl)
end

