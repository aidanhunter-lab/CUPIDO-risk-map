%% Plot satellite data of monthly average phytoplankton quantities from SeaWiFS data

%% Preamble
% Include all required (sub)directories within the search path
project = 'CUPIDO-risk-map';
thisFile = which('map_satellite_phytoplankton_SeaWiFS');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))

%% Map location
coordsTable = readtable('regional bounding coordinates.csv');
disp(coordsTable)

% Choose area from coordsTable -- it's probably useful to vectorise over all areas
area = coordsTable.Location{10};

% Get map coordinates -- these may be specified as input arguments to
% plotBaseMap, and may also be returned as output arguments.
[lon, lat] = plotBaseMap(area, 'createMap', false, 'coordsTable', coordsTable); %, 'coordsTable', coordsTable);


%% Load data
filename = 'SeaWiFS_Phytoplankton-Size-Class-1997-2007-Southern-Ocean.mat';
load(filename) % Stored as variable 'Dat'
disp(Dat)
% Convert -999 to nan
fields = fieldnames(Dat);
for i = 1:length(fields)
    x = Dat.(fields{i});
    x(x == -999) = nan;
    Dat.(fields{i}) = x; clear x
end
% Exclude data outside map boundaries
lae = 0.025 * diff(lat); % allow some data to overlap map boundaries -- avoids white-space
loe = 0.025 * diff(lon);
ind = lat(1) - lae <= Dat.latitude & Dat.latitude <= lat(2) + lae & ...
    lon(1) - loe <= Dat.longitude & Dat.longitude <= lon(2) + lae;

nlon = max(max(sum(ind)));
nlat = max(max(sum(ind,2)));
Dat = structfun(@(z) reshape(z(ind), [nlon, nlat, 12]), Dat, 'UniformOutput', false);


%% Plot maps

% The SeaWiFS global satellite data provides chlorophyll concentration
% estimates fractioned by cell size. Let's use some conversions to change
% from chl units into carbon units...
unitConversions = readtable('Unit Conversions.csv');

convertRatio = 'r_p_cchl'; % plankton carbon/chl
group = 'OpenWater_median';

r_p_cchl = unitConversions.Value(...
    strcmp(unitConversions.Parameter, convertRatio) & strcmp(unitConversions.Group, group));

% Assume the same conversion factor for plankton size classes, and estimate
% carbon concentrations.
fields = fieldnames(Dat);
ind = contains(fields, 'Tchl'); % index the chlorophyll measurements
for i = 1:length(fields)
    if ind(i)
        field = fields{i};
        x = r_p_cchl .* Dat.(field);
        field = strrep(field, 'Tchl', 'carbon');
        Dat.(field) = x;
    end
end

Var = 'Tchl'; % variable to plot
% Var = 'carbon'; % variable to plot

% SeaWiFS data are monthly means => plot separate map for each month
Months.name = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
Months.nameFull = {'January','February','March','April','May','June','July','August','September','October','November','December'};
Months.index = 1:12;
Months.summerMonths = {'Oct','Nov','Dec','Jan','Feb','Mar'};
Months.winterMonths = {'Apr','May','Jun','Jul','Aug','Sep'};
Months.issummer = ismember(Months.name, Months.summerMonths);
Months.iswinter = ismember(Months.name, Months.winterMonths);

multiPanel = false; % all maps on same figure?
switch multiPanel, case true
    nrows = floor(nmonths ^ 0.5);
    ncols = ceil(nmonths / nrows);
end

% Choose colour map
cm = viridis; % viridis is appropriate for phytoplankton
% other options include: plasma, inferno, magma, fake_parula, parula

% Use same colour scale across different maps => not all colours present on
% maps for different months
nlevel = size(cm, 1) - 1;
% colRescale = @(z) log10(z); % use log-scale to produce nice colour gradient
% invColScale = @(z) 10 .^ z;

% Create a separate figure for each month. Store these in a struct so they
% may be combined into a single multipanel figure. Saving figures for
% individual months may be useful for creating animations -- gif files to
% include in presentations to show spatio-temporal dynamics of prey field.

displayMonth = true;
abbrevMonth = false;
monthLat = max(lat) + 0.05 * diff(lat); % text position
monthLon = -150; % bottom-left
monthSize = 13;
monthAlignHoriz = 'right';
monthAlignVert = 'top';

mainTitle = true;
titleText = {'Phytoplankton concentration: 1997-2007 average', 'SeaWiFS global data'};
titleLat = max(lat) + 0.05 * diff(lat); % title position
titleLon = 0; % centred

titleSize = 13;
fullTitle = true; % false => abbreviated month
cbarTitleSize = 12; % text sizing for colourbar
cbarLabelSize = 9;
axisSize = 9;
landColour = .4 .* ones(1,3);
XaxisLocation = 'top';
nColourBackTicks = 7; % number of ticks on colourbar
plotSize = [6 7]; % a little extra height to accomodate the multi-line title

anyData = nan(12,1);

% Make map plot for each month
for m = 1:12
    month = Months.index(m);
    if m == 1, clearvars plotHandles plotAxes plotColourAxes bc; end
    plotName = ['plot' num2str(m)];
    plotHandles.(plotName) = figure;
    set(plotHandles.(plotName), {'Units', 'Position'}, {'inches', [0 0 plotSize(1) plotSize(2)]})
    % Extract data for month m
    ind = Dat.month == month;
    dat = structfun(@(z) reshape(z(ind), [size(z, 1), size(z, 2)]), Dat, 'UniformOutput', false);
    z = dat.(Var);
    % Create map
    plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
        'redrawCoastline', false, 'XaxisLocation', XaxisLocation, 'axesLabelSize', axisSize);
    hold on
    anyData(m) = any(~isnan(z(:)));
%     if anyData
        m_pcolor(dat.longitude, dat.latitude, z)
        set(gca, 'ColorScale', 'log')
        % m_pcolor(dat.longitude, dat.latitude, z_)
%     end
    clear dat    
    plotAxes.(plotName) = gca;
    set(plotAxes.(plotName), 'Colormap', cm)
    
    % redraw map grid
    plotBaseMap(area, 'coordsTable', coordsTable, 'edgecolour', landColour, ...
        'redrawCoastline', false, 'XaxisLocation', XaxisLocation, ...
        'Xticklabels', [], 'axesLabelSize', axisSize);

    % Create colourbar -- full colour range for each figure, sort out scales later...
    if ~anyData(m)
        plotColourAxes(m,:) = nan(1,2);
    else
        plotColourAxes(m,:) = caxis;  % store colour bar limits for each axis
    end
    cb.(plotName) = colorbar(plotAxes.(plotName), 'SouthOutside');
    set(cb.(plotName), 'FontSize', cbarLabelSize);
    switch Var
        case 'Tchl'
            set(cb.(plotName).Label, {'String', 'FontSize'}, ...
                {'total chlorophyll (mg / m^3)', cbarTitleSize});
        case 'carbon'
            set(cb.(plotName).Label, {'String', 'FontSize'}, ...
                {'total carbon (mg / m^3)', cbarTitleSize});
    end

    switch displayMonth, case true
        switch abbrevMonth
            case true, pm = Months.name(m);
            case false, pm = Months.nameFull(m);
        end
        m_text(monthLon, monthLat, pm, ...
            'FontSize', monthSize, ...
            'HorizontalAlignment', monthAlignHoriz, 'VerticalAlignment', monthAlignVert)
    end

    switch mainTitle, case true
        m_text(titleLon, titleLat, titleText, ...
            'FontSize', titleSize, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom')
    end

    pause(0.5)
end

% Rescale colours so that colour bars are identical between monthly maps.
% This should reduce the colour range displayed for some monthas that do
% not attain very high, or low, values.
disp(plotColourAxes)
CLim = [0.1, 10]; % values < CLim(1) are darkest shade; values > CLim are lightest shade
for i = 1:12
    gc = plotAxes.(['plot' num2str(i)]);
    set(gc, 'CLim', CLim)
    cbi = cb.(['plot' num2str(i)]);
    ticks = cbi.TickLabels;
    ticks{1} = ['<' ticks{1}];
    ticks{end} = ['>' ticks{end}];
    set(cbi, 'TickLabels', ticks)
end

% save plots
for m = 1:12
    mon = num2str(Months.index(m));
    if length(mon) == 1, mon = ['0' mon]; end
    filename = ['SeaWiFS_' area '_' Var '_' mon '.png'];
%     filename = ['SeaWiFS_' area '_' Var '_' Months.name{m} '.png'];
    filepath = fullfile(baseDirectory, 'MatLab', 'plots', filename);
    plt = plotHandles.(['plot' num2str(m)]);
    exportgraphics(plt, filepath)
end


GOT TO HERE SO FAR, WHICH IS FINE FOR INDIVIDUAL PLOTS TO DISPLAY ON A POSTER...
    THE CODE BELOW IS OLD AND STILL NEEDS ADJUSTED.









% For South Georgia, all months excluding June & July have data. However,
% there appear to be some spurious measures in May, so let's also exclude
% May.
plotHandles.plot5 = []; plotHandles.plot6 = []; plotHandles.plot7 = [];
plotAxes.plot6 = []; plotAxes.plot7 = [];
cb.plot6 = []; cb.plot7 = [];
plotColourAxes(5:7,:) = nan;

% Create 2 or 5 column multipanel plot to display all the data.
% But first, match the colour scales between months...

% plotHandles
% plotAxes
% plotColourAxes
% cb

pca = 10 .^ plotColourAxes;

switch Var
    case 'carbon'
        caxx = colRescale(...
            [ceil(max(pca(:,1), [], 'omitnan') / 10) * 10, floor(mean(pca(:,2), 'omitnan') / 100) * 100]);
    case 'Tchl'
        caxx = colRescale(...
            [max(pca(:,1), [], 'omitnan'), mean(pca(:,2), 'omitnan')]);
end

cp = 0.05 * diff(caxx);
minz_ = cp + caxx(1); maxz_ = -cp + caxx(2);
% Set tick marks
md = 2; % min displayed digits per tick label
x_ = linspace(minz_, maxz_, 15);
x = 10.^ x_;
y = cellstr(string(x));
y = cellfun(@(z) z(1:strfind(z, '.') - 1), y, 'UniformOutput', false);
ndigits = cellfun(@(z) length(z), y);
y = round(x, 1, 'significant');
ndigits(x < 1) = cell2mat(cellfun(@(z) length(num2str(z)), num2cell(y(x < 1)), 'Uniformoutput', false));
rdigits = ndigits;
rdigits(y < 1) = 2;
for r = 1:length(x), x(r) = round(x(r), rdigits(r), 'significant'); end
switch Var
    case 'carbon'
        x = round(x / 5) * 5; % round to nearest 5
    case 'Tchl'
        x = round(x * 20) / 20; % round to nearest 0.5
end

x = x([1,5,8,12,15]);
Ticks = log10(x);
TickLabels = num2cell(x');
TickLabels = cellfun(@(z) num2str(z), TickLabels, 'UniformOutput', false);
TickLabels{1} = ['<' TickLabels{1}];
TickLabels{end} = ['>' TickLabels{end}];

for m = 1:12
    plotName = ['plot' num2str(m)];
    if isempty(plotAxes.(plotName)), continue; end
    caxis(plotAxes.(plotName), caxx)
    set(cb.(plotName), {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
end

% Save separate plots for each month
savePlots = false;
dirSlides = fullfile(dirBase, 'meetings', 'ecosystem group', 'Science Open Day 29_03_2022', 'animation');
numericTag = true; % false => file name contains abbreviated month; true => numeric tag
nlabs = 1:12;
nlabs = nlabs([7:12, 1:6]); % change numerical order to start midyear/winter
nlabs = cellstr(string(nlabs));
for m = 1:12, if length(nlabs{m}) == 1, nlabs{m} = ['0' nlabs{m}]; end; end

switch savePlots, case true
    for m = 1:nmonths
        plotName = ['plot' num2str(m)];
        if isempty(plotHandles.(plotName)), continue; end
        switch Var
            case 'Tchl', filename = ['SeaWiFS_total_Chl_' strrep(area,' ','_') '_'];
            case 'carbon', filename = ['SeaWiFS_total_carbon_' strrep(area,' ','_') '_'];
        end
        switch numericTag
            case false, filename = [filename Months.name{m} '.png'];
            case true, filename = [filename nlabs{m} '.png'];
        end
%         filepath = fullfile(dirData, filename);
        filepath = fullfile(dirSlides, filename);
        exportgraphics(plotHandles.(plotName), filepath)
    end
end


% All data on sinlge multipanel plot
nrows = 3;
ncols = 3;
jj = [7:12 1:6];
switch multiPanel, case true
    if isfield(plotHandles, 'multiPanelPlot'), plotHandles = rmfield(plotHandles, 'multiPanelPlot'); end
    j = 0;
    for m = jj
        plotName = ['plot' num2str(m)];
        if isempty(plotHandles.(plotName))
            continue
        else
            j = j +1;
        end
        if ~isfield(plotHandles, 'multiPanelPlot')
            plotHandles.multiPanelPlot = figure;
        end
        axCopy = copyobj(plotAxes.(plotName), plotHandles.multiPanelPlot);
        s(j) = subplot(nrows, ncols, j, axCopy);
    end
end
mcb = colorbar(s(8), 'SouthOutside');
set(mcb, {'Ticks', 'TickLabels', 'FontSize'}, {Ticks, TickLabels, cbarLabelSize})
s(8).Position(2) = s(9).Position(2);
s(8).Position(4) = s(9).Position(4);
mcb.Position(1) = s(7).Position(1) + 0.5 * s(7).Position(3);
mcb.Position(3) = (s(9).Position(1) + 0.5 * s(9).Position(3)) - (s(7).Position(1) + 0.5 * s(7).Position(3));
switch Var
    case 'Tchl'
        set(mcb.Label, {'String', 'FontSize'}, ...
            {'total chlorophyll (mg / m^3)', cbarTitleSize});
    case 'carbon'
        set(mcb.Label, {'String', 'FontSize'}, ...
            {'total carbon (mg / m^3)', cbarTitleSize});
end
suptitle(titleMultiPanel)
set(plotHandles.multiPanelPlot, {'Units', 'Position'}, {'inches', [0 0 10 8]})

% adjust panel sizes to reduce white space -- this didn't work well!
exh = 0.0; % extra horizontal space outwith  plot bounds
ew = 1.0; % width expansion
for m = 1:length(s)
    cc = 1 + mod(m - 1, ncols); % column number
    xw = ew * (1 + 2 * exh) / ncols; % panel widths
    xw_ = repmat(xw, [1, ncols]);
    xp = linspace(-exh, 1 + exh - xw, ncols); % bottom left corners
    s(m).Position(1) = xp(cc);
    s(m).Position(3) = xw_(cc);
end

% save multipanel plot

% Save map plot
switch Var
    case 'Tchl'
        filename = ['SeaWiFS_total_Chl_' area '.png'];
    case 'carbon'
        filename = ['SeaWiFS_total_carbon_' area '.png'];
end
filepath = fullfile(dirData, filename);
exportgraphics(plotHandles.multiPanelPlot, filepath)

close all




% % Highlight areas on map? This will be useful for highlighting regions in
% % coordsTable on a map of the entire Southern Ocean, perhaps then comparing
% % chl between regions.
% highlightRegions = false; % Skip the following lines if false
% switch highlightRegions
%     case true
%         regionCoords = readtable('Southern Ocean region coordinates.csv');
%         highlightRegionNames = {'Western Antarctic Peninsula', 'Scotia Sea'};
%         % highlightRegionNames = {'Western Antarctic Peninsula', 'Scotia Sea', 'Bransfield Strait', 'South Georgia'};
%         highlightRegionColour = [1, 0, 0];
%         nregions = length(highlightRegionNames);
%         % Store highlight areas as lists of lon-lat coords of bounding-box vertices
%         coordVertices = cell(1, nregions);
%         for j = 1:nregions
%             region = highlightRegionNames{j};
%             % Adjust coordinate bounding-box for Antarctic Peninsula so that the box
%             % captures more coastline...
%             switch region
%                 case 'Western Antarctic Peninsula'
%                     x = [-75; -65; -65; -60.5; -60; -58; -58; -65; -75; -75];
%                     y = [-70; -70; -66; -64; -64; -63.3; -61; -61; -64; -70];
%                 otherwise
%                     d = regionCoords(strcmp(regionCoords.Location, region),:);
%                     x = [d.Longitude_min; d.Longitude_min; d.Longitude_max; d.Longitude_max; d.Longitude_min];
%                     y = [d.Latitude_min; d.Latitude_max; d.Latitude_max; d.Latitude_min; d.Latitude_min];
%             end
%             coordVertices{:,j} = [x, y];
%         end
% end





close all

% Compare mean chlorophyll concentrations between months, within selected
% regions
plotVals = cell(nregions, nmonths);
for r = 1:nregions
    d = coordVertices{:,r};
    for m = 1:nmonths
        month = months(m);
        ind = Dat.month == month;
        dat = structfun(@(z) reshape(z(ind), [size(z, 1), size(z, 2)]), Dat, 'UniformOutput', false);
        %     d = regionCoords(strcmp(regionCoords.Location, region),:);
        ind = inpolygon(dat.longitude, dat.latitude, d(:,1), d(:,2));
        dat = structfun(@(z) z(ind), dat, 'UniformOutput', false);
        plotVals{r,m} = dat.(Var);
    end
end
nvals = cellfun(@(z) length(z), plotVals(:,1));
mvals = max(nvals);
for j = 1:nregions
    if nvals(j) < mvals
        p = plotVals(j,:);
        n = mvals - nvals(j);
        for jj = 1:length(p)
            p{jj} = [p{jj}; nan(n, 1)];
        end
        plotVals(j,:) = p;
    end
end


plotVals = plotVals(:)';
plotVals = cell2mat(plotVals);

pltRegions = repmat(highlightRegionNames', [nmonths, 1]);
pltMonths = reshape(repmat(Months.name, [nregions, 1]), [nmonths * nregions, 1]);

G1 = repmat(pltMonths, [1, size(plotVals, 1)])';
G2 = repmat(pltRegions, [1, size(plotVals, 1)])';
G1 = categorical(G1(:), Months.name);

pltMeans = figure;
boxchart(G1, plotVals(:), 'GroupByColor', G2(:))
set(gca, 'YScale', 'log')
legend('Location', 'south')
ylabel('Mean chlorophyll concentration (mg / m^3)')

filename = 'SeaWiFS_total Chl distributions_ScotiaSea and AntarcticPeninsula.png';
filepath = fullfile(dirData, filename);
exportgraphics(pltMeans, filepath)


%% Examine Krillbase data in relation to satellite phytoplankton data
dirKrillBase = fullfile(fileparts(fileparts(which('map_satellite_phytoplankton_SeaWiFS'))), ...
    'data', 'KRILLBASE', 'density data');
addpath(dirKrillBase)

filename = 'krillbase_data.csv';
krill = readtable(filename);

r = 'Western Antarctic Peninsula';
cv = coordVertices{1,strcmp(highlightRegionNames, r)};
krill_peninsula = krill(inpolygon(krill.LONGITUDE, krill.LATITUDE, cv(:,1), cv(:,2)),:);
krill_peninsula.region = repmat({r}, [height(krill_peninsula), 1]);
r = 'Scotia Sea';
cv = coordVertices{1,strcmp(highlightRegionNames, r)};
krill_ScotiaSea = krill(inpolygon(krill.LONGITUDE, krill.LATITUDE, cv(:,1), cv(:,2)),:);
krill_ScotiaSea.region = repmat({r}, [height(krill_ScotiaSea), 1]);
krill_ = [krill_peninsula; krill_ScotiaSea];

% Include years/months present in satellite data averages
[yr, mo, day] = datevec(krill_.DATE);
krill_.year = yr;
krill_.monthIndex = mo;
krill_.day = day;
ind = ismember(krill_.monthIndex, Months.index);
krill_ = krill_(ind,:);
[yr, mo, day] = datevec(krill_.DATE);
mi = repmat(Months.index, [height(krill_), 1]);
mo_ = repmat(mo, [1, length(Months.index)]);
x = mo_ == mi;
[I, J] = find(x);
[~,oI] = sort(I);
krill_.month = Months.name(J(oI))';
years = 1997:2007;
nyears = length(years);
ind = ismember(krill_.year, years);
krill_ = krill_(ind,:);
% Use day- or night-time samples
dayornight = 'day';
ind = strcmp(krill_.DAY_NIGHT, dayornight);
krill_ = krill_(ind,:);
% Omit zeros
krill_ = krill_(krill_.STANDARDISED_KRILL_UNDER_1M2 > 0,:);

% Plot krill density data for each month
cols = cm(round(linspace(1, size(cm, 1), nyears)),:); % colour points by sample year
pltKrillDensity = figure;
set(pltKrillDensity, {'Units', 'Position'}, {'inches', [0 0 12 6]})
for r = 1:nregions
    subplot(1, nregions, r)
    region = highlightRegionNames{r};
    dat = krill_(strcmp(krill_.region, region),:);
    % Group density measures by month & year
    plotVals = cell(nyears,nmonths);
    plotCols = cell(nyears,nmonths);
    for y = 1:nyears
        year = years(y);
        indy = dat.year == year;
        if any(indy)
            for m = 1:nmonths
                month = Months.name{m};
                indm = strcmp(dat.month, month);
                if any(indm)
                    d = dat(indy & indm,:);
                    plotVals{y,m} = d.STANDARDISED_KRILL_UNDER_1M2;
                end
            end
        end
    end
    pv = plotVals(:);
    pm = max(cellfun(@(z) length(z), pv));
    for j = 1:length(pv), pv{j} = [pv{j}; nan(pm - length(pv{j}), 1)]; end
    plotVals = reshape(pv, size(plotVals));
    plotVals = plotVals(:)';
    plotVals = cell2mat(plotVals); % columns ordered as years/months (1997Oct,1998Oct,...)
    pltYears = repmat(years', [nmonths, 1]);
    pltMonths = reshape(repmat(Months.name, [nyears, 1]), [nmonths * nyears, 1]);
    G1 = repmat(pltMonths, [1, size(plotVals, 1)])';
    G2 = repmat(pltYears, [1, size(plotVals, 1)])';
    G1 = categorical(G1(:), Months.name);
    bc = boxchart(G1, plotVals(:), 'GroupByColor', G2(:));
    for j = 1:length(bc)
        bc(j).BoxFaceColor = cols(j,:);
        bc(j).MarkerColor = cols(j,:);
    end
    set(gca, 'YScale', 'log')
    if r == 1, legend('Location', 'southwest'); end
    if r == 1, ylabel('krill abundance (individuals / m^2)'); end
    title(region)
    if r == nregions, suptitle('Krill density measurements from KRILLBASE'); end
end


% Match the KrillBase measurements to chlorophyll measures by finding the
% closest points in terms of lon-lats.
% The chlorophyll data are monthly averages over a ten year period, whereas
% the krill density measurements are single samples.
% However, there may still be correlation between chlorophyll and krill
% densities, especially if the satellite chlorophyll data are sensitive to
% persistant (over years) hotspots.
krill_.Tchl = nan(height(krill_),1);
ck = [krill_.LONGITUDE, krill_.LATITUDE, krill_.monthIndex];
mkrill = unique(ck(:,3));
for j = 1:length(mkrill)
    ind = Dat.month == mkrill(j);
    dat = structfun(@(z) reshape(z(ind), [size(ind, 1), size(ind, 2)]), Dat, 'UniformOutput', false);
    indk = ck(:,3) == mkrill(j);
%     kdat = krill_(indk,:);
    nk = sum(indk);
    lon2d = repmat(dat.longitude(:,1), [1, nk]);
    lat2d = repmat(dat.latitude(1,:)', [1, nk]);
    x = abs(lon2d - reshape(ck(indk,1), [1, nk]));
    x = x == min(x);
    [I,~] = find(x);
    x = abs(lat2d - reshape(ck(indk,2), [1, nk]));
    x = x == min(x);
    [J,~] = find(x);
    z = dat.Tchl(I,J);
    krill_.Tchl(indk) = diag(z);
end

f1 = figure; set(f1, {'Units', 'Position'}, {'inches', [0 0 6 6]})
scatter(krill_.Tchl, krill_.STANDARDISED_KRILL_UNDER_1M2)
xlabel('Mean chlorophyll concentration (mg / m^3)')
ylabel('krill abundance (individuals / m^2)')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

f2 = figure; set(f2, {'Units', 'Position'}, {'inches', [0 0 6 6]})
scatter(krill_.CLIMATOLOGICAL_TEMPERATURE, krill_.Tchl)
xlabel(['Krill samples climatological temperature (', char(176), 'C)'])
ylabel('Mean chlorophyll (mg / m^3)')
% set(gca, 'YScale', 'log')

f3 = figure; set(f3, {'Units', 'Position'}, {'inches', [0 0 6 6]})
scatter(krill_.CLIMATOLOGICAL_TEMPERATURE, krill_.STANDARDISED_KRILL_UNDER_1M2)
xlabel(['Krill samples climatological temperature (', char(176), 'C)'])
ylabel('krill abundance (individuals / m^2)')
set(gca, 'YScale', 'log')




% The remaining variables are relative abundances => can all be expessed in
% the same units (proportions: 0 <= p <= 1) => the maps should/could share
% colour scales.
Vars = {'Micro_percent_Tchl', 'Nano_percent_Tchl', 'Pico_percent_Tchl'};
z = [];
for i = 1:length(Vars)    
    z_ = Dat.(Vars{i});
    z = [z; z_(:)];
end
nlevel = size(cm, 1);
zl = quantile(z(:), linspace(0, 1, nlevel+1)); % colour scale levels

multiPanel = false; % all maps on same figure?

cm = plasma; % colour map
% colormap viridis
% colormap plasma
% colormap inferno
% colormap magma
% colormap fake_parula
% colormap parula

% Make two map types: (1) a separate multipanel plot for each month where
% each panel displays a size class; (2) a separate multipanel plot for each
% size class where each panel displays a month.


for m = 1:length(months)
    month = months(m);
    ind = Dat.month == month;
    
    dat = structfun(@(z) reshape(z(ind), [size(z, 1), size(z, 2)]), Dat, 'UniformOutput', false);
    
    
    
    
    
    
    
    
    z = dat.(Var);
    
    mapName = ['month' num2str(m)];
    
    switch multiPanel
        case true
            if m == 1
                map = figure;
                set(map, {'Units', 'Position'}, {'inches', [0 0 16 12]})
            end
            subplot(nrows, ncols, m)
        case false
            maps.(mapName) = figure;
            set(maps.(mapName), {'Units', 'Position'}, {'inches', [0 0 8 8]})
    end
    
    plotBaseMap(area, 'edgecolour', .4 .* ones(1,3), 'redrawCoastline', false);
    hold on
    if any(~isnan(z(:)))
        [~,mc] = m_contourf(dat.longitude, dat.latitude, z, zl);
        set(mc, 'LineStyle', 'none')
    end
    
    mapAxes.(mapName) = gca;
    set(mapAxes.(mapName), 'Colormap', cm)
    pause(0.5)
end
























%% Map entire continent
area = 'continent';
mapContinent = figure;
set(mapContinent, {'Units', 'Position'}, {'inches', [0 0 8 8]})
[lon, lat, rad] = plotBaseMap_(area, 'edgecolour', .4 .* ones(1,3), 'redrawCoastline', false);
% The output variables to plotBaseMap are optional -- a convenience for
% filtering data later...
% plotBaseMap(area, 'edgecolour', .4 .* ones(1,3), 'redrawCoastline', false)
% plotBaseMap(area, 'rad', 40, 'areacolour', [0.15, 0.75, 0.1], 'edgecolour', [0.25, 0.5, 0.2], 'redrawCoastline', true)


%% Map peninsula
area = 'peninsula';
mapPeninsula = figure;
set(mapPeninsula, {'Units', 'Position'}, {'inches', [0 0 8 8]})
% axes('position', [0, 0, 1, 0.8]) % leave space for legend
[lon, lat] = plotBaseMap_(area, 'edgecolour', .4 .* ones(1,3), 'redrawCoastline', false);
% plotBaseMap(area, 'edgecolour', .4 .* ones(1,3), 'redrawCoastline', false)
% plotBaseMap(area, 'lon', [-75,-25], 'lat', [-70, -50], 'redrawCoastline', true)

%% Load data
filename = 'SeaWiFS_Phytoplankton-Size-Class-1997-2007-Southern-Ocean.mat';
load(filename)
disp(Dat)
% Convert -999 to nan
fields = fieldnames(Dat);
for i = 1:length(fields)
    x = Dat.(fields{i});
    x(x == -999) = nan;
    Dat.(fields{i}) = x; clear x
end
% Exclude data outside map boundaries
switch area
    case 'continent'
        ind = Dat.latitude >= lat & Dat.latitude <= lat + rad;
    case 'peninsula'
        ind = lat(1) <= Dat.latitude & Dat.latitude <= lat(2) & ...
            lon(1) <= Dat.longitude & Dat.longitude <= lon(2);
end
measures = {'Tchl','Micro_percent_Tchl','Nano_percent_Tchl','Pico_percent_Tchl'};
for i = 1:length(measures)
    x = Dat.(measures{i});
    x(~ind) = nan;
    Dat.(measures{i}) = x;
end

m = 1;
ind = Dat.month == 1;

dat = structfun(@(z) reshape(z(ind), [size(z, 1), size(z, 2)]), Dat, 'UniformOutput', false);

z = dat.Tchl;
z = log10(z);
z = z - min(z(:));
tmp = figure;
histogram(z(:))
nlevel = 11;
zl(2:nlevel+1) = [linspace(quantile(z(:), 0.05), quantile(z(:), 0.95), nlevel - 1), ...
    max(z(:))];
close(tmp)
hold on
[~,mc] = m_contourf(dat.longitude, dat.latitude, z, zl);
% [~,mc] = m_contourf(dat.longitude, dat.latitude, log10(dat.Tchl), zl);
set(mc, 'LineStyle', 'none')












% Plot points
hold on
pointSize = 40;

% LEGEND DESIGN MUST BE DIFFERENT FOR PLOTTING MULTIPLE YEARS OF KRILL DATA
% datType
% ai
n = length(pltData.Latitude);

switch ai
    case 'onlyKrill'
        x = [pltData.Year, strcmp(pltData.DayOrNight, 'night')];
        [~, I] = sortrows(x);
        pd = pltData;
        pltData = rmfield(pltData, 'dat');
        pltData = structfun(@(z) z(I,:), pltData, 'UniformOutput', false);
        pltData.dat = pd.dat; clear pd
        x = [pltData.Year, strcmp(pltData.DayOrNight, 'day')];
        x = [zeros(1,2); diff(x)];
        legi = (x(:,1) | x(:,2));
    case 'multiple'
        t1 = table(pltData.Source);
        t2 = table(unique(pltData.Source)); t2.index = (1:height(t2))';
        x = join(t1, t2);
        x = [0; diff(x.index)];
        legi = x ~= 0;
end

clear ppleg
for i = 1:n
    pp = m_scatter(pltData.Longitude(i), pltData.Latitude(i));    
    pp.MarkerEdgeColor = pltData.Colour(i,:);
    pp.SizeData = pointSize;
    pp.Marker = pltData.Shape{i};
    % create dummy points for legend
    if i == 1
        ppleg(1) = pp;
    end
    if legi(i)
        ppleg(length(ppleg)+1) = pp;
    end
end

% Legend
legPosition = 'eastoutside';
legOrientation = 'vertical';
leg = legend(ppleg, 'Location', legPosition, 'Orientation', legOrientation);

switch ai
    case 'multiple'
        x = v.datSource;
        t = v.datType;
        g = v.datGroup;        
        for i = 1:length(x)
            y = x{i};
            if contains(y, '_')
                y = strsplit(y, '_');
                y = [y{1} ' (' y{2} ')'];
            end
            x{i} = y;
            if ~strcmp(y, 'density data')
                x{i} = [x{i} ': ' t{i}];
            else
                x{i} = [t{i} ': ' g{i}];
            end
        end
        set(leg, 'String', x)
    case 'onlyKrill'
        legi(1) = true;
        y = pltData.Year(legi);
        dn = pltData.DayOrNight(legi);
        g = v.datGroup{1};
        clear x
        for i = 1:sum(legi)
            x{i} = [num2str(y(i)) ': ' dn{i} '-time'];
        end
        set(leg, 'String', x);
        set(leg.Title, 'String', [g ': abundance'])
end


% Save map plot
switch area
    case 'continent', map = mapContinent;
    case 'peninsula', map = mapPeninsula;
end
switch ai
    case 'multiple'
        filename = ['map_Antarctic ' area '_sampleLocations.png'];
        folder = fullfile(fileparts(which('map_Antarctic')), 'plots', 'maps');
        filepath = fullfile(folder, filename);
        exportgraphics(map, filepath)
    case 'onlyKrill'
        filename = ['map_Antarctic ' area '_sampleLocations_KRILLBASE_2010-2016.png'];
        folder = fullfile(fileparts(which('map_Antarctic')), 'plots', 'maps');
        filepath = fullfile(folder, filename);
        exportgraphics(mapPeninsula, filepath)
end



