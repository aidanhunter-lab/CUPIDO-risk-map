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
thisFile = which('download_chl.m');
dirBase = fileparts(fileparts(thisFile));
dirMatLab = fullfile(dirBase, 'MatLab');
dirData = fullfile(dirBase, 'data', 'chlorophyll', 'SeaWiFS');
addpth(dirMatLab)
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

