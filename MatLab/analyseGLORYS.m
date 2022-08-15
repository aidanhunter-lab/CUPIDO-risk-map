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

% Loading data from a single year requires >4GB RAM, so we are limited in
% how much to load into RAM. Let's limit it to a single year or season,
% i.e, 12 months.

% Choose a SINGLE season from 1993-2019. (A season is defined as starting 
% July 1st and ending June 30th.)
season = 2019;
nseasons = length(season);
ymin = season;
ymax = season + 1;
yrs = ymin:ymax;
% Choose month(s), within seasons
% months = [7:12, 1:6]; % Jul-Dec, Jan-Jun
months = [10:12, 1:3];
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
clearvars imonth month iyr yr season nseasons ymin ymax yrs months nmonths
disp(timeDat)

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
    clearvars f filename fileInfo dataDim
end

% Convert the time variable to MatLab standard time
GLORYS_startDate = [1950 1 1 0 0 0];
adjustTime = datenum(GLORYS_startDate);
dat.time = dat.time / 24; % days since 1/1/1950
dat.time = dat.time + adjustTime;
% datevec(double(dat.time))
clearvars GLORYS_startDate adjustTime i


%% Summary plots

savePlots = true;

% Water density will vary most strongly with depth. It may also vary
% significantly with latitude and, probably to a lesser extent, with
% longitude.

% Visualise how density varies with depth
months = timeDat.month(1);
imonths = ismember(timeDat.month, months);

v = permute(dat.density, [3 4 1 2]); % plotting variable -- permute depth to first dimension
v = v(:,imonths,:,:);
mv = mean(v, [2, 3, 4], 'omitnan'); % Average over all lat-lons and months
sv = std(v, 0, [2, 3, 4], 'omitnan'); % Find variability across all lat-lons and months
er = [mv - sv, mv + sv];

figure
plot(mv, dat.depth, 'Color', [0 0 0])
hold on
plot(er, dat.depth, ':', 'Color', [0 0 0])
set(gca, 'YDir', 'reverse')
xlabel('Seawater density (kg m^{-3})')
ylabel('Depth (m)')
% Seawater density varies very little with respect to latitude, longitude 
% or month -- depth is by far the main source of variability.



% Map density against lat-lon, for single depth layers.
coordsTable = readtable('regional bounding coordinates.csv');
disp(coordsTable)
area = coordsTable.Location{10};
[~, lat] = plotBaseMap(area, 'createMap', false, 'coordsTable', coordsTable);

% Choose month(s) and depth(s) to plot
months = timeDat.month(1:end);
nmonths = length(months);
imonths = ismember(timeDat.month, months);

% The modelled depths form several discrete layers, so if I choose an
% arbitrary depth then the model output should be interpolated to match
% that depth.
depth_range = [min(dat.depth), max(dat.depth)];
disp(depth_range) % choose some values within this range
% Given that water density varies much more strongly with depth than with
% other variables choose a single depth value to plot, otherwise we cannot
% discern variability within individual maps.
depths = 100;
ndepths = length(depths);
multipleDepths = ndepths > 1;
% slice the model output at the chosen month(s)
z = dat.density(:,:,:,imonths);
% interpolate to the chosen depth(s)
zi = nan(nlon, nlat, ndepths, nmonths);
for i = 1:nmonths
    zi(:,:,:,i) = interp3(dat.lat, dat.lon, dat.depth, z(:,:,:,i), dat.lat, dat.lon, single(depths));
end


% Plot
landColour = [0.4 0.4 0.4];
XaxisLocation = 'top';
axisSize = 9;
cbarLabelSize = 9;
cbarTitleSize = 12; % text sizing for colourbar
plotSize = [5 7]; % a little extra height to accomodate the multi-line title
% plotSize = [6 6]; % a little extra height to accomodate the multi-line title
mainTitle = true;
% titleText = {'GLORYS model output', 'Seawater density'}; % multiline title
titleText = 'GLORYS model output'; % multiline title
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred
titleSize = 13;

displayDepth = true; % show the depth/month
displayMonth = true; % as row/column labels
monthLabAdj = 0.15;

switch multipleDepths
    case false
        npanels = nmonths * ndepths;
        ncols = ceil(sqrt(npanels));
        nrows = ceil(npanels / ncols);
        rows = repmat(1:nrows, [ncols, 1]);
        cols = repmat(reshape(1:ncols, [ncols, 1]), [1, nrows]);
        multiPanel = npanels > 1;
    case true
        % If multiple depths or months are selected then create a multipanel plot
        % with depth along y-direction and month along x-direction
        nrows = ndepths;
        ncols = nmonths;
        rows = repmat(1:nrows, [ncols, 1]);
        cols = repmat(reshape(1:ncols, [ncols, 1]), [1, nrows]);
        npanels = nrows * ncols;
        multiPanel = npanels > 1;
end

cm = inferno; % Choose colourmap (see the utility directory for nice options)

switch multiPanel
    case true
        plt = figure;
        set(plt, {'Units', 'Position'}, {'inches', [0 0 ncols * plotSize(1) nrows * plotSize(2)]})
%         plt.Clipping = 'off';
        plotColourAxes = nan(npanels,2);
        for i = 1:npanels
            row = rows(i);
            col = cols(i);
            subplot(nrows,ncols,i)
            axisName = ['axis', num2str(i)];
            plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
                'redrawCoastline', false, 'XaxisLocation', XaxisLocation, 'axesLabelSize', axisSize);
            hold on
            plotAxes.(axisName) = gca;
            set(plotAxes.(axisName), 'Colormap', cm)
            pos_orig = plotAxes.(axisName).Position;
            switch multipleDepths
                case false
                    z_ = zi(:,:,:,i); % data slice to plot
                case true
                    z_ = zi(:,:,row,col);
            end
            z_ = permute(z_, [2, 1]);
            m_pcolor(dat.lon, dat.lat, z_)
            % redraw map grid
            plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
                'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ...
                'Xticklabels', [], 'axesLabelSize', axisSize);
            % Labels
            switch displayMonth, case true
                m_abbrev = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
                m_name = {'January','February','March','April','May','June','July','August','September','October','November','December'};
                switch multipleDepths
                    case true
                        if row == 1
                            Title = m_abbrev{months(col)};
                            m_text(0, max(lat) + monthLabAdj * diff(lat), Title, ...
                                'HorizontalAlignment', 'center')
                            %                     title(Title, 'FontWeight', 'normal')
                        end
                    case false
                        Title = m_abbrev{months(i)};
                        m_text(0, max(lat) + monthLabAdj * diff(lat), Title, ...
                            'HorizontalAlignment', 'center')
                end
                 pos_orig(4) = 0.9 * pos_orig(4); % make a little extra room for titles
            end
            switch displayDepth, case true
                % Only makes sense when multiple depths are plotted
                switch multipleDepths, case true
                    if col == 1
                        d_lab = [num2str(depths(row)) ' m deep'];
                        ylabel(d_lab)
                    end
                end
            end

           
            % Colour bars
            cb.(axisName) = colorbar(plotAxes.(axisName), 'SouthOutside', ...
                'Visible', 'off');
            plotColourAxes(i,:) = caxis; % store colour bar limits for each axis
            % Colourbar moves the map. Reset in original position
            plotAxes.(axisName).Position = pos_orig;

            pause(0.25)

        end

        % Rescale colours so that colour bars are identical between plot panels.
        % This will reduce the colour range displayed for some panels that do
        % not attain very high, or low, values.

        cp = [0.1, 0.9];
        cq = nan(2,ndepths);
        for i = 1:ndepths
            cq(:,i) = quantile(reshape(zi(:,:,i,:), 1, []), cp);
        end
        CLim = [min(cq(1,:)), max(cq(2,:))];
        for i = 1:npanels
            axisName = ['axis' num2str(i)];
            set(plotAxes.(axisName), 'CLim', CLim)
        end

        % Create a single colour bar showing scale for all plot panels
        colbar = cb.axis1;
        evenColumns = mod(ncols, 2) == 0; % positioning depends on number of columns
        switch evenColumns
            case true
                p1 = (nrows-1) * ncols + 0.5 * ncols;
                p2 = p1 + 1;
                pos1 = cb.(['axis' num2str(p1)]).Position;
                pos2 = cb.(['axis' num2str(p2)]).Position;
                colbar.Position(1) = pos1(1);
                colbar.Position(3) = pos2(1)+pos2(3)-pos1(1);
                colbar.Position(2) = pos1(2);
                colbar.Position(4) = pos1(4);
            case false
                if ncols > 1
                    p1 = (nrows-1) * ncols + 0.5 * (ncols+1) - 1;
                    p2 = p1 + 1; p3 = p1 + 2;
                    pos1 = cb.(['axis' num2str(p1)]).Position;
                    pos2 = cb.(['axis' num2str(p2)]).Position;
                    pos3 = cb.(['axis' num2str(p3)]).Position;
                    colbar.Position(1) = 0.5 * ((pos1(1) + 0.5 * pos1(3)) + (pos2(1) + 0.5 * pos2(3)));
                    colbar.Position(3) = 0.5 * ((pos2(1) + 0.5 * pos2(3)) + (pos3(1) + 0.5 * pos3(3))) - colbar.Position(1);
                    colbar.Position(2) = pos1(2);
                    colbar.Position(4) = pos1(4);
                elseif ncols == 1
                    colbar.Position = cb.(['axis' num2str(nrows)]).Position;
                end
        end
        colbar.Visible = 'on';

        set(colbar, 'FontSize', cbarLabelSize);
        set(colbar.Label, {'String', 'FontSize'}, ...
            {'Seawater density (kg m^{-3})', cbarTitleSize});
        
        switch mainTitle, case true
            switch multipleDepths
                case false
                    sgtitle([titleText ': ' num2str(depths) ' m deep']);
                case true
                    sgtitle(titleText)
            end
        end
end



switch savePlots, case true
    filename = 'GLORYS seawater density maps.png';
    filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
    exportgraphics(plt, filepath)
end




