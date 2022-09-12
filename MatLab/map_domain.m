%% Specify model domain and resolution
function v = map_domain(varargin)

% Optional arguments include:
% Longitude & latitude ranges has 2-element vectors: lon_range & lat_range
% Resolution in longitudinal & latitudianl directions: dlon & dlat.
% Calling the function without optional arguments returns default values.
extractVarargin(varargin)

%% Horizontal domain

% Domain boundaries
if ~exist('lon_range', 'var')
    lon_range = [-180, 180];
end
if ~exist('lat_range', 'var')
    lat_range = [-90, -50];
end

% Horizontal resolution
if ~exist('dlon', 'var')
    dlon = 9; % degrees longitude
end
if ~exist('dlat', 'var')
    dlat = 3; % degrees latitude
end

% Test for gaps in grid -- the resolution should fit neatly into the range
lon_diff = diff(lon_range);
lat_diff = diff(lat_range);
lon_mod = mod(lon_diff, dlon);
lat_mod = mod(lat_diff, dlat);
lon_fit = lon_mod == 0;
lat_fit = lat_mod == 0;
% Adjust domain for gaps...
% Shrink latitude range to match resolution
if ~lat_fit
    lat_range(2) = lat_range(2) - lat_mod;
end
% Alter longitude resolution to match range
if ~lon_fit
    dlon = dlon:-1:1;
    x = mod(lon_diff, dlon) == 0;
    dlon = dlon(find(x, 1));
end

% Store domain in output struct
v.lon_range = lon_range;
v.lat_range = lat_range;
v.dlon = dlon;
v.dlat = dlat;
v.nlon = diff(lon_range) / dlon; % number of longitudinal
v.nlat = diff(lat_range) / dlat; % and latitudinal grid cells
v.longrid = (lon_range(1):dlon:lon_range(2))'; % grid cell boundaries
v.latgrid = (lat_range(1):dlat:lat_range(2))';
v.lon = 0.5 * (v.longrid(1:end-1) + v.longrid(2:end)); % grid cell midpoints
v.lat = 0.5 * (v.latgrid(1:end-1) + v.latgrid(2:end));


%% Vertical domain
% Most microbial remineralsing activity occurs above 500m. Dynamics below
% 1000m are less interesting, although we do want to know whether faecal
% pellets reach the seafloor. Thus, model the upper water at higher
% resolution than lower water column, and have a coarse resolution below
% 1000m while keeping seafloor depth within the model.

if ~exist('maxDepth', 'var')
    maxDepth = 6000; % maximum modelled (seafloor) depth
end
depth_range = [0, maxDepth];

if ~exist('microAttDepth', 'var')
    microAttDepth = 500; % depth where microbial action attenuates to become negligible
end

if ~exist('exportDepth', 'var')
    exportDepth = 1000; % depth at which sinking matter is considered exported to depth
end

if ~exist('surfaceLayerWidth', 'var')
    surfaceLayerWidth = 20; % width of uppermost modelled depth layer
end

if ~exist('nDepthsUpper', 'var')
    nDepthsUpper = 8; % number of depth layers in upper water column, above microAttDepth
end

if ~exist('deepLayerWidth', 'var')
    deepLayerWidth = 1000; % width of modelled depth layers below the export depth
end

constraints = false(1,2);
constraints(1) = nDepthsUpper <= microAttDepth / surfaceLayerWidth; % required so that upper water column depth layers have widths increasing with depth

depthLayerGradient = (microAttDepth - nDepthsUpper*surfaceLayerWidth) / ... 
    sum(1:nDepthsUpper-1);
depthLayerWidths = surfaceLayerWidth + depthLayerGradient .* (0:nDepthsUpper-1);

d = exportDepth - microAttDepth;
constraints(2) = d >= 0;
if d > 0
    depthLayerWidths(nDepthsUpper+1) = d;
end


maxDepth = ceil(maxDepth / deepLayerWidth) * deepLayerWidth; % round up to accomodate choice of depth layer width
nDeepLayers = (maxDepth - exportDepth) / deepLayerWidth; % maximum number of modelled depth layers below the export depth

depthLayerWidths(nDepthsUpper+2:nDepthsUpper+nDeepLayers+1) = repmat(deepLayerWidth, [1, nDeepLayers]);
depthEdges = cumsum([0, depthLayerWidths]);
depthMid = 0.5 * (depthEdges(1:end-1) + depthEdges(2:end));

v.depth_range = depth_range;
v.depthgrid = depthEdges';
v.depth = depthMid';


% Area at top and bottom of each grid cell
longrid = repmat(v.longrid, [1, v.nlat+1]);
latgrid = repmat(reshape(v.latgrid, 1, []), [v.nlon+1, 1]);
Area = areaquad(latgrid(1:end-1,1:end-1), longrid(1:end-1,1:end-1), ...
    latgrid(2:end,2:end), longrid(2:end,2:end));
h = 4 * pi .* (earthRadius - v.depthgrid) .^ 2;
h = repmat(reshape(h, 1, 1, []), v.nlon, v.nlat, 1);
AreaEarth = h .* Area;
v.area = AreaEarth;

% Volume of each grid cell
Area = repmat(Area, [1 1 length(v.depth)]);
Volume = 4 / 3 * pi .* Area;
h = abs(diff((earthRadius - v.depthgrid) .^ 3));
VolumeEarth = Volume .* repmat(reshape(h, 1, 1, []), v.nlon, v.nlat, 1);
v.volume = VolumeEarth;


if ~constraints(1)
    warning('Depth layers in upper water column have widths that decrease with depth! It should be the other way around. Reduce values of inputs nDepthsUpper or surfaceLayerWidth.')
end

if ~constraints(2)
    warning('Constraint exportDepth >= microAttDepth is not satisfied -- alter these input variables. This constraint may be removed in future...')
end
