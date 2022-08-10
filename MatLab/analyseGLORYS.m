% Analyse GLORYS physical model -- summary statistics, make plots

%% Preamble

% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('analyseGLORYS.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))
dataPath = fullfile(baseDirectory, 'data', 'physical_models', ... 
    'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern Ocean');
Months = ["01","02","03","04","05","06","07","08","09","10","11","12"];
Days = ["15","16"];

%% Load data

% Choose a year (range)
ymin = 2019;
ymax = 2019;
yrs = ymin:ymax;
nyrs = length(yrs);
% Choose month(s)
% months = [12 1 2];
months = 1:12;
nmonths = length(months);
ntimes = nyrs * nmonths;

% Try loading the water density data for the Southern Ocean for an entire
% year.

for i = 1:nyrs
    for j = 1:nmonths
        fileCount = (i - 1) * nyrs + j;
        f = {fullfile(dataPath, append(num2str(yrs(i)), Months(months(j) == str2double(Months)), Days(1), ".nc")), ...
            fullfile(dataPath, append(num2str(yrs(i)), Months(months(j) == str2double(Months)), Days(2), ".nc"))}; % possible file names for this year & month
        filename = f{[exist(f{1}, 'file'), exist(f{2}, 'file')] == 2}; % choose the file name that exists in the search path
        if i == 1 && j == 1
%             ncdisp(filename)
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
        dat.time(fileCount) = ncread(filename, 'time');
        dat.density(:,:,:,fileCount) = ncread(filename, 'density');
    end
end

% Convert the time variable to MatLab standard time
GLORYS_startDate = [1950 1 1 0 0 0];
adjustTime = datenum(GLORYS_startDate);
dat.time = dat.time / 24; % days since 1/1/1950
dat.time = dat.time + adjustTime;
% datevec(double(dat.time))


%% Summary plots

% Map density against lat-lon, for single depth layers or averaged over depth range

coordsTable = readtable('regional bounding coordinates.csv');
disp(coordsTable)
area = 'Southern Ocean';
icoords = strcmp(coordsTable.Location, area);
lon = [coordsTable.Longitude_min(icoords), coordsTable.Longitude_max(icoords)];
lat = [coordsTable.Latitude_min(icoords), coordsTable.Latitude_max(icoords)];

% Choose month and depth to plot
month = 3;
% The modelled depths form several discrete layers, so if I choose an
% arbitrary depth then the model output should be interpolated to match
% that depth... try interp3 (X, Y, Z = lon, lat, depth)

landColour = [0.4 0.4 0.4];
XaxisLocation = 'top';
axisSize = 9;
cbarLabelSize = 9;
cbarTitleSize = 12; % text sizing for colourbar
plotSize = [6 7]; % a little extra height to accomodate the multi-line title
mainTitle = true;
titleText = {'GLORYS model output', 'Seawater density'}; % multiline title
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred
titleSize = 13;

plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, 'axesLabelSize', axisSize);
hold on
plotAxes = gca;
pos_orig = plotAxes.Position;
z2D = permute(dat.density(:,:,1,1), [2 1]);
m_pcolor(dat.lon, dat.lat, z2D)
% redraw map grid
plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
    'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ...
    'Xticklabels', [], 'axesLabelSize', axisSize);
CLim = quantile(z2D(:), [0.2, 0.8]);
set(gca, 'CLim', CLim)
colormap inferno
cb = colorbar(gca, 'SouthOutside');
set(cb, 'FontSize', cbarLabelSize);
set(cb.Label, {'String', 'FontSize'}, ...
    {'Seawater density (kg m^{-3})', cbarTitleSize});
% Colourbar moves the map. Reset in original position
plotAxes.Position = pos_orig;

% Set colourbar width to 90% of map width
cbWidth = 0.9 * (pos_orig(3) - pos_orig(1));
cb.Position(1) = 0.5 - 0.5 * cbWidth;
cb.Position(3) = cbWidth;
switch mainTitle, case true
    m_text(titleLon, titleLat, titleText, ...
        'FontSize', titleSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom')
end






% Density vs depth (and latitude)

md = mean(dat.density, [1 2]); 
scatter(squeeze(md(:,:,:,1)), dat.depth)

md = mean(dat.density); % mean over longitude

help contourfps

figure
plt = contourfps(repmat(reshape(dat.lat, [1 nlat]), [nlon 1]), repmat(dat.lon, [1 nlat]), dat.density(:,:,1,1), linspace(1024, 1028, 5));
colorbarps





