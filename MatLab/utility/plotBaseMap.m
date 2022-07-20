function varargout = plotBaseMap(area, varargin)
% Returns figure displaying map of study area.
% The study area should be one of those listed in the coordsTable optional
% input argument (see file 'regional bounding coordinates.csv')

extractVarargin(varargin)

if ~exist('areacolour', 'var'), areacolour = [.7, .7, .7]; end
if ~exist('edgecolour', 'var'), edgecolour = [1, 1, 1]; end
if ~exist('XaxisLocation', 'var'), XaxisLocation = 'top'; end
if ~exist('redrawCoastline', 'var'), redrawCoastline = false; end
if ~exist('createMap', 'var'), createMap = true; end
if ~exist('projection', 'var'), projection = 'lambert'; end % default map projection compatible with rectangular coordinate bounds
if ~exist('lon', 'var'), lon = []; end
if ~exist('lat', 'var'), lat = []; end
if ~exist('origin', 'var'), origin = []; end
if ~exist('parallels', 'var'), parallels = []; end
if ~exist('ellipsoid', 'var'), Ellipsoid = []; else, Ellipsoid = ellipsoid; end
if ~exist('axesLabelSize', 'var'), axesLabelSize = 11; end

if strcmp(area, 'Southern Ocean') || strcmp(area, 'Antarctic Continent')
    projection = 'Stereographic';
end


if ~exist('coordsTable', 'var')
    % It's better to pass coordsTable as a variable because the saved table
    % is more likely to be up-to-date than this...
    Location = {'Lazarev Sea'; 'Cosmonaut Sea'; 'Cooperation Sea, Prydz Bay'; 'Mawson Sea'; 'Somov Sea'; 'Western Antarctic Peninsula'; 'Bransfield Strait'; 'Scotia Sea'; 'South Georgia'};
    Longitude_min = [-5    30    70    95   160   -75   -61   -50   -40]';
    Longitude_max =  [5    40    80   105   170   -65   -58   -40   -34]';
    Latitude_min = [-70  -70  -70  -66  -71  -70  -64  -61  -55.5]';
    Latitude_max = [-57   -57   -57   -57   -58   -57   -62   -53   -53.5]';
    coordsTable = table(Location, Longitude_min, Longitude_max, Latitude_min, Latitude_max);
end

areaIndex = ismember(coordsTable.Location, area); % Is the specified area matched to some coordinates?
knownArea = any(areaIndex);

switch knownArea
    case false
        warning('Specified area does not correspond to known coordinates. It is best to run this fuction with optional input argument coordsTable which stores area names and coordinates, otherwise see coordinates hard-coded into plotBaeMap.m')
        return
    case true
        nout = nargout;
        varargout = cell(1,nout);
        Coords = coordsTable(areaIndex,:);
        if isempty(lon), lon = [Coords.Longitude_min, Coords.Longitude_max]; end
        if isempty(lat), lat = [Coords.Latitude_min, Coords.Latitude_max]; end
        if nout > 0
            coords = [lon; lat];
            for i = 1:min(nout,2), varargout{i} = coords(i,:); end
        end
        if isempty(origin)
            origin = [Coords.Originx, Coords.Originy];
            if any(isnan(origin)), origin = []; end
        end
        if isempty(parallels)
            parallels = [Coords.Parallel1, Coords.Parallel2];
            if any(isnan(parallels)), parallels = []; end
        end
        if isempty(Ellipsoid)            
            Ellipsoid = Coords.Ellipsoid{1};
            if isnan(Ellipsoid), Ellipsoid = []; end
        end
        
        oi = [isempty(origin) && isempty(parallels) && isempty(Ellipsoid), ...
            ~isempty(origin) && isempty(parallels) && isempty(Ellipsoid), ...
            isempty(origin) && ~isempty(parallels) && isempty(Ellipsoid), ...
            isempty(origin) && isempty(parallels) && ~isempty(Ellipsoid), ...
            ~isempty(origin) && ~isempty(parallels) && isempty(Ellipsoid), ...
            ~isempty(origin) && isempty(parallels) && ~isempty(Ellipsoid), ...
            isempty(origin) && ~isempty(parallels) && ~isempty(Ellipsoid), ...
            ~isempty(origin) && ~isempty(parallels) && ~isempty(Ellipsoid)];
        orientation = {'NoOriginNoParallelsNoEllipsoid', 'OriginNoParallelsNoEllipsoid', ...
            'NoOriginParallelsNoEllipsoid', 'NoOriginNoParallelsEllipsoid', ...
            'OriginParallelsNoEllipsoid', 'OriginNoParallelsEllipsoid', ...
            'NoOriginParallelsEllipsoid', 'OriginParallelsEllipsoid'};
        orientation = orientation{oi};
end

switch createMap
    case true
        area_ = strrep(area, ' ', '');
        area_ = strrep(area_,',','');
        mapDataFile = ['MapFile_' area_ '_' projection 'Projection.mat'];
        mapFileExists = exist(mapDataFile, 'file') == 2;
        if ~mapFileExists || redrawCoastline
            createMapFiles(area, projection, ...
                'lon', lon, 'lat', lat, ...
                'orientation', orientation, ...
                'parallels', parallels, 'origin', origin, 'Ellipsoid', Ellipsoid, ...
                'areacolour', areacolour, 'edgecolour', edgecolour)
        end
        
        if ~strcmp(projection, 'Stereographic')
            switch orientation
                case 'NoOriginNoParallelsNoEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat);
                case 'OriginNoParallelsNoEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'origin', origin);
                case 'NoOriginParallelsNoEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'parallels', parallels);
                case 'NoOriginNoParallelsEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'ellipsoid', Ellipsoid);
                case 'OriginParallelsNoEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'origin', origin, 'parallels', parallels);
                case 'OriginNoParallelsEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'origin', origin, 'ellipsoid', Ellipsoid);
                case 'NoOriginParallelsEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'parallels', parallels, 'ellipsoid', Ellipsoid);
                case 'OriginParallelsEllipsoid'
                    m_proj(projection, 'lon', lon, 'lat', lat, 'origin', origin, 'parallels', parallels, 'ellipsoid', Ellipsoid);
            end
        else
            m_proj(projection, 'lon', mean(lon), 'lat', -90, 'rad', diff(lat))
        end

        % plot coastline
        m_usercoast(mapDataFile, 'patch', areacolour, 'EdgeColor', edgecolour);
        
        % plot grid
        if ~exist('Xticklabels', 'var') && ~exist('Yticklabels', 'var')
            m_grid('XaxisLocation', XaxisLocation, 'fontsize', axesLabelSize);
        elseif exist('Xticklabels', 'var') && exist('Yticklabels', 'var')
            xtl = Xticklabels; ytl = Yticklabels;
            m_grid('XaxisLocation', XaxisLocation, ...
                'xticklabels', xtl, 'yticklabels', ytl, 'fontsize', axesLabelSize);
        elseif exist('Xticklabels', 'var') && ~exist('Yticklabels', 'var')
            xtl = Xticklabels;
            m_grid('XaxisLocation', XaxisLocation, 'xticklabels', xtl, 'fontsize', axesLabelSize);
        elseif ~exist('Xticklabels', 'var') && exist('Yticklabels', 'var')
            ytl = Yticklabels;
            m_grid('XaxisLocation', XaxisLocation, 'yticklabels', ytl, 'fontsize', axesLabelSize);
        end
end

