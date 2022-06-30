%% Load and plot CTD data for South Georgian waters
% These data were downloaded from BODC

%% Preamble

% Include all subdirectories within search path
addpath(genpath(fileparts(which('map_CTD_SouthGeorgia'))))

% Set data source and region   Choose region and season(s) to match directories where data are saved
source = 'BODC';
region = 'South Georgia';

%% Load & clean the data

% Map all CTD samples and return data filtered to include only depth profiles

[Data, MetaData, cruises] = cleanData_CTD(source, region, ... 
    'plotMaps', true, 'displayData', false, ...
    'flagGoodMeasure', [49 50 56], 'plotSpuriousPoints', false, ... 
    'removeSpuriousPoints', true, 'onlyDepthProfiles', true, ...
    'ShowMapTitle', true);
% Other output objects could be included in above function...

disp(head(Data))
disp(head(MetaData))

%% Plot the data...

plotVars = {'Chlorophyll', 'Temperature', 'SalinityPractical'}; % must match Data names
plotVars_ = plotVars; plotVars_([1,3]) = {[plotVars_{1} ' a'], 'Salinity'}; % adjust for plot axes labels
unitVars = {'mg/m^3', '\circC', 'unitless'};
nVars = length(plotVars);
colVars = {[0 1 0], [1 0 0], [0 0 1]};

% Separate subplots for each CTD cast -- separate figures for each cruise

ucruises = unique(Data.Cruise, 'stable');
ncruises = length(ucruises);
Rescale = true;

for i = 1:ncruises
    cruiseID = ucruises{i};
    indi = strcmp(Data.Cruise, cruiseID);
    usamples = unique(Data.Label(indi), 'stable');
    nsamples = length(usamples);
    nrows = ceil(nsamples ^ 0.5); ncols = ceil(nsamples / nrows);
    plotHandles.(cruiseID) = figure;
    set(plotHandles.(cruiseID), {'Units', 'Position'}, {'inches', [0, 0, 4*ncols, 4*nrows]})
    for j = 1:nsamples
        subplot(nrows, ncols, j)
        sampleID = usamples{j};
        indj = indi & strcmp(Data.Label, sampleID);
        Dat = Data(indj,:);
        x = Dat{:,plotVars};
        switch Rescale, case true
            x = (x - min(x)) ./ (max(x) - min(x));
        end
        y = -Dat.Depth;
        for k = 1:length(plotVars)
            plot(x(:,k), y, 'Color', colVars{k})
            if k == 1, hold on; end
        end
        title(sampleID, 'FontWeight', 'normal'); ylabel('depth (m)')
    end
    sgtitle(['Depth profiles from cruise ' cruiseID])
%     suptitle(['Depth profiles from cruise ' cruiseID])
end

fields = fieldnames(plotHandles);
for i = 1:length(fields), close(plotHandles.(fields{i})); end; clearvars plotHandles


% Separate figure for each station, including all cruises and variables
ustations = unique(Data.StationID); % find all stations
substations = strfind(ustations, '.'); % and put them in 'correct' sequential order
tt = cellfun(@(z) isempty(z), substations);
ustations_ = cell(size(ustations));
ustations_(tt) = ustations(tt); ustations_ = cellfun(@(z) str2double(z), ustations_);
ustations_ = cellstr(string(ustations_));
for i = 1:sum(tt)
    tt_ = find(tt);
    if tt_(i) >= length(tt), break; end
    j = tt_(i)+1:tt_(i+1)-1;
    usub = ustations(j);
    jj = cellstr(string(sort(str2double(cellfun(@(z) z(strfind(z, '.')+1:end), usub, 'UniformOutput', false)))));
    for j = 1:length(jj)
        jj{j} = [ustations_{tt_(i)} '.' jj{j}];
    end
    ustations_(tt_(i)+1:tt_(i+1)-1) = jj;
end
ustations = ustations_; clearvars ustations_
nstations = length(ustations);

lineTypes = {'-','--',':','-.'}; % indicates different cruises sampling the same station

for i = 1:nstations
    stationID = ustations{i};
    sLab = ['station' strrep(stationID, '.', '_')];
    indi = strcmp(Data.StationID, stationID);
    ucruises = unique(Data.Cruise(indi), 'stable');
    ncruises = length(ucruises);
    ncols = nVars; nrows = 1;
    if nVars > 4
        ncols = ceil(nVars .^ 0.5); nrows = ceil(ncruises / ncols);
    end
    plotHandles.(sLab) = figure;
    set(plotHandles.(sLab), {'Units', 'Position'}, {'inches', [0, 0, 4*ncols, 4*nrows]})
    for j = 1:nVars
        subplot(nrows, ncols, j)
        %         Dat = Data(indi, {'Cruise', 'Depth', plotVars{j}});
        for k = 1:ncruises
            cruiseID = ucruises{k};
            indk = indi & strcmp(Data.Cruise, cruiseID);
            usamples = unique(Data.Label(indk), 'stable');
            nsamples = length(usamples);
            for l = 1:nsamples
                sampleID = usamples{l};
                indl = indk & strcmp(Data.Label, sampleID);
                Dat = Data(indl,:);
                plot(Dat.(plotVars{j}), Dat.Depth, 'Color', colVars{j}, 'LineStyle', lineTypes{k})
%                 plot(Dat.(plotVars{j}), -Dat.Depth, 'Color', colVars{j})
                if k == 1 && l == 1, hold on; end
            end
        end
        set(gca, 'YDir', 'reverse')
        xlabel([plotVars_{j} ' (' unitVars{j} ')']); ylabel('Depth (m)')
    end
    sgtitle(['Station ' stationID])
end


% Create plots that are useful for slides
stationID = '5'; % pick a station -- let's just use station 5
axisSize = 11;
fontSize = 16;
fsm = fontSize / axisSize;
lwd = 2;

sLab = ['station' strrep(stationID, '.', '_')];
indi = strcmp(Data.StationID, stationID);
ucruises = unique(Data.Cruise(indi), 'stable');
ncruises = length(ucruises);
ncols = nVars; nrows = 1;
if nVars > 4
    ncols = ceil(nVars .^ 0.5); nrows = ceil(ncruises / ncols);
end
plt = figure;
set(plt, {'Units', 'Position'}, {'inches', [0, 0, 4*ncols, 4*nrows]})
for j = 1:nVars
    subplot(nrows, ncols, j)
    %         Dat = Data(indi, {'Cruise', 'Depth', plotVars{j}});
    for k = 1:ncruises
        cruiseID = ucruises{k};
        indk = indi & strcmp(Data.Cruise, cruiseID);
        usamples = unique(Data.Label(indk), 'stable');
        nsamples = length(usamples);
        for l = 1:nsamples
            sampleID = usamples{l};
            indl = indk & strcmp(Data.Label, sampleID);
            Dat = Data(indl,:);
            plot(Dat.(plotVars{j}), Dat.Depth, 'Color', colVars{j}, 'LineStyle', lineTypes{k}, 'LineWidth', lwd)
            %                 plot(Dat.(plotVars{j}), -Dat.Depth, 'Color', colVars{j})
            if k == 1 && l == 1, hold on; end
        end
    end
    set(gca, {'YDir', 'FontSize', 'LabelFontSizeMultiplier'}, {'reverse', axisSize, fsm})
    xlabel([plotVars_{j} ' (' unitVars{j} ')'])
    if j == 1
        ylabel('Depth (m)')
        gc = gca; xl = gc.XLim; yl = gc.YLim;
        text(xl(1)-0.0*diff(xl), yl(1) - 0.05 * diff(yl), ...
            ['station ' stationID], 'FontSize', mean([axisSize fontSize]), ...
            'HorizontalAlignment', 'left')
%         text(xl(2)-0.05*diff(xl), yl(2) - 0.05 * diff(yl), ...
%             ['station ' stationID], 'FontSize', mean([axisSize fontSize]), ...
%             'HorizontalAlignment', 'right')
    end
end
% sgtitle(['Station ' stationID])

savePlots = true;
dirBase = fileparts(fileparts(which('map_CTD_SouthGeorgia')));
dirSlides = fullfile(dirBase, 'meetings', 'ecosystem group', 'Science Open Day 29_03_2022');
switch savePlots, case true
    plotName = 'depth profile plot for slides';
    filename = ['CTD_station' stationID '_' plotName '.png'];
    filepath = fullfile(dirSlides, filename);
    exportgraphics(plt, filepath)
end





% Create plots that are useful for slides -- as above but use a single
% panel that may be animated through the different variables

stationID = '5'; % pick a station -- let's just use station 5

sLab = ['station' strrep(stationID, '.', '_')];
indi = strcmp(Data.StationID, stationID);
ucruises = unique(Data.Cruise(indi), 'stable');
ncruises = length(ucruises);
ncols = 1; nrows = 1;
for j = 1:nVars
    plotHandles.(plotVars{j}) = figure;
    set(plotHandles.(plotVars{j}), {'Units', 'Position'}, {'inches', [0, 0, 4*ncols, 4*nrows]})
    subplot('Position', [0.14, 0.12, 0.76, 0.8])
    for k = 1:ncruises
        cruiseID = ucruises{k};
        indk = indi & strcmp(Data.Cruise, cruiseID);
        usamples = unique(Data.Label(indk), 'stable');
        nsamples = length(usamples);
        for l = 1:nsamples
            sampleID = usamples{l};
            indl = indk & strcmp(Data.Label, sampleID);
            Dat = Data(indl,:);
            plot(Dat.(plotVars{j}), Dat.Depth, 'Color', colVars{j})
            %                 plot(Dat.(plotVars{j}), -Dat.Depth, 'Color', colVars{j})
            if k == 1 && l == 1, hold on; end
        end
    end
    axesHandles.(plotVars{j}) = findall(plotHandles.(plotVars{j}),'type','axes');
    set(axesHandles.(plotVars{j}), {'YDir'}, {'reverse'})
%     set(axesHandles.(plotVars{j}), {'YDir', 'Units', 'OuterPosition', 'Position'}, {'reverse', 'inches', [0 0 4 4], [0.75 0.75 3 3]})
    xlabel([plotVars_{j} ' (' unitVars{j} ')']); ylabel('Depth (m)')
end
% for j = 1:nVars
%     if strcmp(plotVars{j}, 'Chlorophyll'), continue; end
% %     set(axesHandles.(plotVars{j}), 'TightInset', axesHandles.Chlorophyll.TightInset)
%     set(axesHandles.(plotVars{j}), 'Position', axesHandles.(plotVars{1}).Position)
% end

% Save individual plots of each variable
savePlots = true;
dirBase = fileparts(fileparts(which('map_CTD_SouthGeorgia')));
dirSlides = fullfile(dirBase, 'meetings', 'ecosystem group', 'Science Open Day 29_03_2022', 'animation2');
switch savePlots, case true
    for i = 1:nVars
        plotName = plotVars{i};
        if isempty(plotHandles.(plotName)), continue; end
        filename = ['CTD_station' stationID '_' plotName '.png'];
        filepath = fullfile(dirSlides, filename);
        exportgraphics(plotHandles.(plotName), filepath)
    end
end



    
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





