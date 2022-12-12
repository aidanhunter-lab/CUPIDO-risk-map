function [out, sources, nsources] = filter_plastic_conc_data2(dat, DAT, varargin)
% Filters plastic abundance data sets and combines into a single table.
% Inputs 'dat' and DAT are a struct and table of data sets produced from
% load_plastic_conc_data.m

%% Optional arguments
extractVarargin(varargin)

if ~exist('onScreenDisplay', 'var')
    onScreenDisplay = false;
end
if ~exist('saveData', 'var')
    saveData = false;
end
if ~exist('outputFilename', 'var')
    filename = 'plastic_quantity.csv';
else
    filename = outputFilename;
end

if ~exist('FilterByLitterCategory', 'var')
    FilterByLitterCategory = false;
end

if ~exist('FilterByLitterScale', 'var')
    FilterByLitterScale = false;
end

if ~exist('regularisePlasticForm', 'var')
    regularisePlasticForm = false;
end

% if ~exist('returnOnlyMicroplasticMeasures', 'var')
%     returnOnlyMicroplasticMeasures = false;
% end

%% Directories
project = 'CUPIDO-risk-map';
thisFile = which('filter_plastic_conc_data.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);

%% Data sources
sources = fieldnames(dat);
make_source_ref = @(DAT) strrep(strrep(strrep(strrep(DAT.Source, ' ', '_'), '(', ''), ')', ''), '-', '_');
source_ref = make_source_ref(DAT);
source_ind = @(source_ref, source) strcmp(source_ref, source); % index DAT by source
DATbySource = @(DAT, source_ref, source) DAT(source_ind(source_ref, source),:); % filter DAT by source
nsources = length(sources);

%% Measurement units used by each source?
for i = 1:nsources
    source = sources{i};
    d = DATbySource(DAT, source_ref, source);
    units.(source) = unique(d.Unit, 'stable')';
end
switch onScreenDisplay, case true, disp(units); end

% Measurement units may depend on sample type (sediment samples may be
% different from water samples).
uSampleTypes = [];
for i = 1:nsources
    source = sources{i};
    d = DATbySource(DAT, source_ref, source);
    uSampleTypes = unique([uSampleTypes; d.SampleType], 'stable');
    sampleTypes.(source) = unique(d.SampleType, 'stable')';
end

% Volumetric concentration 'pieces/m^3' is most common unit for seawater 
% samples; some data use area density 'pieces/km^2', some use mass density
% 'g/km^2'.
varOpts_water = {'concentration', 'massConcentration', 'density', 'massDensity'}; %, 'all'}; % filtering options
varOpts_sediment = {'concentration', 'density'};

Units.water.concentration = {'pieces/m3', 'items/m3', 'particles/m3', 'number/m3'};
Units.water.massConcentration = {'g/m3'};
Units.water.density = {'pieces/km2', 'items/km2', 'number/km2'};
Units.water.massDensity = {'g/km2'};

Units.sediment.concentration = {'pieces/g'};
Units.sediment.density = {'pieces/m2'};

Units_water = struct2cell(Units.water);
useUnit_water = cellfun(@(z) z(1), Units_water);
Units_sediment = struct2cell(Units.sediment);
useUnit_sediment = cellfun(@(z) z(1), Units_sediment);

% Standardise Variable column and omit rows with with units not listed in
% Units -- unwanted Variable
for i = 1:nsources
    source = sources{i};
    sampleType = sampleTypes.(source);
    if contains(sampleType, 'water')
        Units_ = Units_water;
        useUnit = useUnit_water;
        varOpts = varOpts_water;
    end
    if contains(sampleType, 'sediment')
        Units_ = Units_sediment;
        useUnit = useUnit_sediment;
        varOpts = varOpts_sediment;
    end
%     d = DATbySource(source);
    I = find(source_ind(source_ref, source));
    n = length(I);
    for j = n:-1:1
        k = cellfun(@(z) ismember(DAT.Unit{I(j)}, z), Units_);
        if any(k)
            DAT.Variable(I(j)) = varOpts(k);
            DAT.Unit(I(j)) = useUnit(k);
        else
            DAT.Variable(I(j)) = {''};
        end
    end
    DAT(strcmp(DAT.Variable, ''),:) = []; % omit unwanted rows
    source_ref = make_source_ref(DAT);
    dat.(source).abundance = DATbySource(DAT, source_ref, source);
end


%% Regularise the data

for i = 1:nsources
    source = sources{i};
    d = DATbySource(DAT, source_ref, source);
    
    % SampleType
    d.SampleType = lower(d.SampleType);

    % SampleGear
    d.SampleGear = strrep(d.SampleGear, 'Sample bottle', 'bottle');
    d.SampleGear = strrep(d.SampleGear, 'Visual survey', 'visual survey');
    j = contains(d.SampleGear, ' towed neuston nets');
    d.SampleGear = strrep(d.SampleGear, ' towed neuston nets', '');
    d.SampleGear(j) = append(cellstr(repmat('neuston net', sum(j), 1)), ' ', d.SampleGear(j));
    d.SampleGear(j) = strrep(d.SampleGear(j), '-', ' ');
    d.SampleGear = strrep(d.SampleGear, 'mu m', [char(hex2dec('03BC')) 'm']);
    d.SampleGear = strrep(d.SampleGear, 'mum', [char(hex2dec('03BC')) 'm']);

    % PlasticForm
    d.PlasticForm = lower(d.PlasticForm);
    d.PlasticForm = strrep(d.PlasticForm, 'plastic', 'particle');
    d.PlasticForm = strrep(d.PlasticForm, 'microfibre', 'fibre');
    d.PlasticForm = strrep(d.PlasticForm, 'microplastic', 'particle');
    % The plastic forms could be regularised here, but it's probably useful
    % to retain the original descriptions for the interactive map.

    % Plastic size
    d.PlasticSize = lower(d.PlasticSize);
    d.PlasticSize = strrep(d.PlasticSize, 'mu m', [char(hex2dec('03BC')) 'm']);
    d.PlasticSize = strrep(d.PlasticSize, 'mum', [char(hex2dec('03BC')) 'm']);
    
    % SiteCategory
    d.SampleGear = strrep(d.SampleGear, 'Open ocean', 'open ocean');
    
    % SampleID -- set as numeric counter starting 1
    j = d.SampleAtStation; % some data sources contain sample stations and longer tows, these get separate SampleID lists
    x = d(:,'SampleID');
    % samples at fixed stations
    SampleID = unique(x.SampleID(j), 'stable');
    nx = length(SampleID);
    SampleID_new = (1:nx)';
    ux1 = table(SampleID, SampleID_new);
    % samples from longer tows
    SampleID = unique(x.SampleID(~j), 'stable');
    nx = length(SampleID);
    SampleID_new = (1:nx)';
    ux2 = table(SampleID, SampleID_new);
    ux = [ux1; ux2];
    x = join(x, ux);
    d.SampleID = cellstr(num2str(x.SampleID_new));

    % Regularise Depth
    d.Depth = lower(d.Depth);
    d.Depth = strrep(d.Depth, ' m', 'm');

    % Measurement Statistic
    d.Statistic = lower(d.Statistic);
    d.Statistic = strrep(d.Statistic, 'minimum', 'min');
    d.Statistic = strrep(d.Statistic, 'maximum', 'max');
    d.Statistic = strrep(d.Statistic, 'average', 'mean');
    d.Statistic = strrep(d.Statistic, 'stdev', 'sd');

    DAT(source_ind(source_ref, source),:) = d;
    dat.(source).abundance = d;
end

% Filter data by litter category or plastic scale?
LitterCategories = unique(DAT.LitterCategory);
switch FilterByLitterCategory
    case false
    case LitterCategories
        DAT = DAT(strcmp(DAT.LitterCategory, FilterByLitterCategory),:);
        source_ref = make_source_ref(DAT);
        for i = 1:nsources
            dat.(sources{i}).abundance = DATbySource(DAT, source_ref,sources{i});
        end
    otherwise
        warning('Optional argument FilterByLitterCategory does not match litter categories listed in the data -- no filtering was done.')
end

LitterScales = unique(DAT.LitterScale);
switch FilterByLitterScale
    case false
    case LitterScales
        DAT = DAT(strcmp(DAT.LitterScale, FilterByLitterScale),:);
        source_ref = make_source_ref(DAT);
        for i = 1:nsources
            dat.(sources{i}).abundance = DATbySource(DAT, source_ref,sources{i});
        end
    otherwise
        warning('Optional argument FilterByLitterScale does not match litter scales listed in the data -- no filtering was done.')
end

switch regularisePlasticForm, case true
    % Regularise the naming convention for plastic types
    Fibres = {'line', 'fibre', 'filament', 'line/fibre'};
    isFibre = ismember(DAT.PlasticForm, Fibres);
    Fragments = {'fragment', 'flake', 'granule', 'sphere', 'particle', 'pellet'};
    isFragment = ismember(DAT.PlasticForm, Fragments);
    Films = {'film'};
    isFilm = ismember(DAT.PlasticForm, Films);
    isOther = ~ismember(DAT.PlasticForm, [{'all'}, Fibres, Fragments, Films]) & ...
        ~cellfun(@(z) isempty(z), DAT.PlasticForm);
    DAT.PlasticForm(isFibre) = {'fibre'};
    DAT.PlasticForm(isFragment) = {'fragment'};
    DAT.PlasticForm(isFilm) = {'film'};
    DAT.PlasticForm(isOther) = {'other'};
    for i = 1:nsources
        dat.(sources{i}).abundance = DATbySource(DAT, source_ref,sources{i});
    end
end


% % Some sample dates are not specific, but rather give an interval using 2
% % columns. Therefore all Date variables need 2 columns
% % NO LONGER NEEDED?
% for i = 1:nsources
%     source = sources{i};
%     d = dat.(source).abundance;
%     if size(d.Date, 2) == 1
%         d.Date = [d.Date, NaT(height(d), 1)];
%         dat.(source).abundance.Date = d.Date;
%     end
% end



% % Some data sources do not report the measurement means, but only the min 
% % and max measurements at each station. This is usually because the data
% % were reported as within some range [min, max]. Let's estimate mean values
% % for those data where only min/max is reported simply by taking midpoints.
% 
% % 'M NOT KEEN ON THIS... SHOULD PROBABLY FIND A BETTER WAY TO DISPLAY THE
% % DATA AS IT IS REPORTED BY EACH SOURCE
% 
% % useStat = 'max';
% useStat = 'mean';
% for i = 1:nsources
%     source = sources{i};
%     d = dat.(source).abundance;
%     stats = unique(d.Stat, 'stable');
%     useStatExists = ismember(useStat, stats);
%     if useStatExists
%         d = d(strcmp(d.Stat, useStat),:);
%         dat.(source).abundance = d;
%         continue
%     else
%         switch useStat
%             case 'mean'
%                 j = ismember({'min','max'}, stats);
%                 if all(j)
%                     stations = unique(d.Station);
%                     nstations = length(stations);
%                     inc = height(d) / nstations;
%                     d_ = d(inc:inc:height(d),:);
%                     d_.Stat = repmat({useStat}, height(d_), 1);
%                     for k = 1:length(stations)
%                         x = d.Value(d.Station == k & ismember(d.Stat, {'min','max'}));
%                         if x(1) > 0
%                             x = exp(0.5 * sum(log(x)));
%                         elseif x(1) == 0
%                             x = 0.5 * x(2);
%                         end
%                         d_.Value(d_.Station ==k) = x;
%                     end
%                     d = d_;
%                 else
%                     warning(['Data from ' source ' has been omitted because useStat is unmatched and the Stat column does not contain rows labelled "min" and "max".'])
%                 end
%             case 'max'
%                 j = strcmp(d.Stat, 'mean');
%                 if ~any(j)
%                     warning(['Data from ' source ' has been omitted because useStat is unmatched and the Stat column does not contain rows labelled "mean".'])
%                 end
%                 d = d(j,:);
%         end
%         dat.(source).abundance = d;
%     end
% end % The abundance data should now contain a single measurement per station



% % Omit variables not shared by all data tables -- these occur after the 
% % Value column, and could be included later...
% % THIS NEEDS MORE WORK... GROUPING VARIABLES TO DISTINGUISH PLASTIC
% % CATEGORIES
% for i = 1:nsources
%     source = sources{i};
%     d = dat.(source).abundance;
%     d = d(:,1:find(strcmp(d.Properties.VariableNames, 'Value')));
%     dat.(source).abundance = d;
% end



% for i = 1:nsources
%     source = sources{i};
%     d = dat.(source);
%     s = size(d.abundance, 2);
%     disp([source ' ' num2str(s)])
% end

% % Combine all plastic abundance data into a single table
% out = structfun(@(z) z.abundance, dat, 'UniformOutput', false);
% out = struct2cell(out);
% out = vertcat(out{:});
% 
% sources = unique(out.Source, 'stable');
% nsources = length(sources);

out = DAT;

% Save full table (and separate tables for different measurements
% variables)
switch saveData, case true
    writetable(out, fullfile(baseDirectory, 'data/plastic_quantity', filename))
end

