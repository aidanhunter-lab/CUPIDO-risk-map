function createMapFiles(area, projection, varargin)
%% Create and save map data using GSHHS files

% Call this function from plotBaseMap.m

% These mapping data are good for detailed maps of small areas, but are not
% worthwhile for large areas

extractVarargin(varargin)

if ~exist('orientation', 'var')
    Orientation = '';
else
    Orientation = orientation;
end

area_ = strrep(area, ' ', '');
area_ = strrep(area_, ',', '');

fileName = ['MapFile_' area_ '_' projection 'Projection'];

% create projection
switch projection

    case {'lambert', 'Lambert Conformal Conic'}

        switch Orientation
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
            otherwise
                m_proj(projection, 'lon', lon, 'lat', lat);
        end

    case {'Stereograpic'}
        m_proj(projection, 'lon', mean(lon), 'lat', -90, 'rad', diff(lat))
end


% view map
m_gshhs_h('patch', areacolour, 'edgecolor', edgecolour);
% m_grid();
m_grid('XaxisLocation', 'top');

% save map data file
thisFile = which('createMapFiles.m');

m_gshhs_h('save', fullfile(fileparts(thisFile), fileName));
