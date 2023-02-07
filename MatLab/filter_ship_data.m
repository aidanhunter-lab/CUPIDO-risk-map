%% Filter AIS ship taffic data
% Estimate annual ship time per vessel, total time and time spent within
% specified regions (lon-lat grid).

saveOutput = true;

%% Directories
thisFile = which('filter_ship_data.m');
project = 'CUPIDO-risk-map';
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
dataDirectory = fullfile(baseDirectory, 'data/shipping/McCarthy_2022');

%% Load data
filename = 'so_obs_aidan.csv';
filepath = fullfile(dataDirectory, filename);
data = readtable(filepath);

%% Basic filtering/cleaning -- get data dimensions
dates = datevec(data.date_time);
yrs = unique(dates(:,1));
for i = length(yrs):-1:1
    records_per_year(2,i) = sum(dates(:,1) == yrs(i));
    records_per_year(1,i) = yrs(i);
end
disp(records_per_year)
% There is only a single record for 2013, so omit this.
data(dates(:,1) == 2013,:) = [];
dates = datevec(data.date_time);
yrs = unique(dates(:,1));
nyrs = length(yrs);

all_vessel_id = unique(data.vessel_id, 'stable');
nvessels = length(all_vessel_id);

% Each listed vessel appears in Southern Ocean in one or more years
vessel_observed = false(nvessels, nyrs); % presence/absence matrix indicating, for each year, if vessel was observed
x = unique([data.vessel_id, dates(:,1)], 'rows');
for i = 1:nvessels
    vessel_id = all_vessel_id(i);
    y = x(x(:,1) == vessel_id, 2); % years that vessel i was observed
    vessel_observed(i, ismember(yrs, y)) = true;
end
clearvars x y i

[~, I] = sort(data.move_id);
data = data(I,:);

%% Estimate ship-time

shipTime_total = zeros(nvessels, nyrs);
for i = 1:nvessels
    vessel = all_vessel_id(i);
    ind_i = data.vessel_id == vessel;
    for j = 1:nyrs
        yr = yrs(j);
        wasObserved = vessel_observed(i,j);
        if ~wasObserved, continue; end
        ind_j = ind_i & dates(:,1) == yr;
        dat = data(ind_j,:); % subset data for vessel i, year j
        nd = height(dat);
        if nd == 1, continue; end % A single AIS record is insufficient to determine visit duration
        gaps = diff(dat.move_id) ~= 1; % index when vessels leave the Southern Ocean
        gaps_ind = find(gaps);
        gaps_start = [1; gaps_ind + 1];
        gaps_end = [gaps_ind; nd];
        % Omit visits with a single AIS ping as these data cannot be used 
        % to determine visit duration
        singlePing = gaps_start == gaps_end;
        if any(singlePing)
%             disp(['i=' num2str(i) ': j=' num2str(j)])
            dat(gaps_start(singlePing),:) = [];
            nd = height(dat);
            gaps = diff(dat.move_id) ~= 1; % index when vessels leave the Southern Ocean
            gaps_ind = find(gaps);
            gaps_start = [1; gaps_ind + 1];
            gaps_end = [gaps_ind; nd];
        end
        alltimes = datenum(dat.date_time);
        nvisits = sum(gaps) + 1;
        v = 0;
        while v < nvisits
            v = v + 1;
            visit_ind = false(nd, 1);
            visit_ind(gaps_start(v):gaps_end(v)) = true;
            alltimes_v = alltimes(visit_ind);
            visitTime = alltimes_v(end) - alltimes_v(1);
            shipTime_total(i,j) = shipTime_total(i,j) + visitTime;
        end
    end
end


% Calculate time each ship is within each cell of selected lat-lon grid.
lonres = 9;
latres = 3;
longrid = -180:lonres:180; % grid cell edge positions
latgrid = -90:latres:-60;
nlon = length(longrid) - 1;
nlat = length(latgrid) - 1;
ncells = nlon * nlat;
cell_index = reshape(1:ncells, nlon, nlat);

lonmin = longrid(1:end-1);
lonmax = longrid(2:end);
latmin = latgrid(1:end-1);
latmax = latgrid(2:end);
shipTime_gridded = zeros(nvessels, nyrs, ncells);

for i = 1:nvessels
    vessel = all_vessel_id(i);
    ind_i = data.vessel_id == vessel;
    for j = 1:nyrs
        yr = yrs(j);
        wasObserved = vessel_observed(i,j);
        if ~wasObserved, continue; end
        ind_j = ind_i & dates(:,1) == yr;
        dat = data(ind_j,:); % subset data for vessel i, year j
        nd = height(dat);
        if nd == 1, continue; end % A single AIS record is insufficient to determine visit duration
        gaps = diff(dat.move_id) ~= 1; % index when vessels leave the Southern Ocean
        gaps_ind = find(gaps);
        gaps_start = [1; gaps_ind + 1];
        gaps_end = [gaps_ind; nd];
        % Omit visits with a single AIS ping as these data cannot be used 
        % to determine visit duration
        singlePing = gaps_start == gaps_end;
        if any(singlePing)
            dat(gaps_start(singlePing),:) = [];
            nd = height(dat);
            gaps = diff(dat.move_id) ~= 1; % index when vessels leave the Southern Ocean
            gaps_ind = find(gaps);
            gaps_start = [1; gaps_ind + 1];
            gaps_end = [gaps_ind; nd];
        end
        alltimes = datenum(dat.date_time); % time of each AIS ping
        dlon = dat.longitude;
        dlat = dat.latitude;
        which_lon = lonmin < dlon & dlon <= lonmax;
        which_lat = latmin < dlat & dlat <= latmax;
        which_lon = repmat(which_lon, 1, 1, nlat);
        which_lat = repmat(reshape(which_lat, nd, 1, nlat), 1, nlon, 1);
        cell_index_ = repmat(reshape(cell_index, 1, nlon, nlat), nd, 1, 1);
        where = cell_index_(which_lon & which_lat); % the grid cell occupied during each AIS ping
        nvisits = sum(gaps) + 1;
        moves = diff(where) ~= 0 | gaps; % index when vessels move between grid cells
        moves_ind = find(moves);
        moves_start = [1; moves_ind + 1];
        moves_end = [moves_ind; nd];
        moves = [moves_start moves_end];
        v = 0;
        while v < nvisits
            v = v + 1;
            visit_ind = false(nd, 1);
            visit_ind(gaps_start(v):gaps_end(v)) = true;
            moves_v = moves(moves(:,1) == gaps_start(v) | ...
                moves(:,2) == gaps_end(v),:);
            nplaces = size(moves_v, 1);
            p = 0;
            while p < nplaces
                p = p + 1;
                moves_p = moves_v(p,:);
                place_ind = false(nd, 1);
                place_ind(moves_p(1):moves_p(2)) = true;
                alltimes_p = alltimes(place_ind);
                visitTime = alltimes_p(end) - alltimes_p(1); % time within grid cell
                % Account for times taken to move between grid cells
                if 1 < p && p < nplaces
                    visitTime = visitTime + ...
                        0.5 * (alltimes_p(1) - alltimes(moves_p(1) - 1)) + ...
                        0.5 * (alltimes(moves_p(2) + 1) - alltimes_p(end));
                elseif p == 1 && nplaces > 1
                    visitTime = visitTime + ...
                        0.5 * (alltimes(moves_p(2) + 1) - alltimes_p(end)); % time to move to other cell
                elseif p == nplaces && nplaces > 1
                    visitTime = visitTime + ...
                        0.5 * (alltimes_p(1) - alltimes(moves_p(1) - 1));
                end
                where_ = unique(where(place_ind));
                if length(where_) > 1, error('Error: location indexing problem. We need to index a single grid cell but multiple are selected.'); end
                shipTime_gridded(i,j,where_) = shipTime_gridded(i,j,where_) + visitTime;
            end
        end
    end
end


% Create a table of results to save as a csv file
vessel_id = repmat(all_vessel_id, 1, nyrs, ncells);
vessel_id = vessel_id(:);

vessel_name_u = unique(data(:,["vessel_id", "vessel_name"]), 'stable');
vessel_name_tmp = table(vessel_id);
vessel_name_tmp = join(vessel_name_tmp, vessel_name_u);
vessel_name = vessel_name_tmp.vessel_name;


year = repmat(reshape(yrs, 1, nyrs), nvessels, 1, ncells);
year = year(:);

lonmin_ = repmat(reshape(lonmin, 1, 1, nlon), nvessels, nyrs, 1, nlat);
lonmin = lonmin_(:);
lonmax_ = repmat(reshape(lonmax, 1, 1, nlon), nvessels, nyrs, 1, nlat);
lonmax = lonmax_(:);
latmin_ = repmat(reshape(latmin, 1, 1, 1, nlat), nvessels, nyrs, nlon, 1);
latmin = latmin_(:);
latmax_ = repmat(reshape(latmax, 1, 1, 1, nlat), nvessels, nyrs, nlon, 1);
latmax = latmax_(:);

ship_time = shipTime_gridded(:);

output_9x3 = table(vessel_id, vessel_name, year, lonmin, lonmax, latmin, latmax, ship_time);



% Repeat on a higher resolution grid
lonres = 3;
latres = 1;
longrid = -180:lonres:180; % grid cell edge positions
latgrid = -90:latres:-60;
nlon = length(longrid) - 1;
nlat = length(latgrid) - 1;
longrid2d = repmat(reshape(longrid, nlon + 1, 1), 1, nlat + 1);
latgrid2d = repmat(reshape(latgrid, 1, nlat + 1), nlon + 1, 1);
ncells = nlon * nlat;
cell_index = reshape(1:ncells, nlon, nlat);

lonmin = longrid(1:end-1);
lonmax = longrid(2:end);
latmin = latgrid(1:end-1);
latmax = latgrid(2:end);
shipTime_gridded = zeros(nvessels, nyrs, ncells);

for i = 1:nvessels
    vessel = all_vessel_id(i);
    ind_i = data.vessel_id == vessel;
    for j = 1:nyrs
        yr = yrs(j);
        wasObserved = vessel_observed(i,j);
        if ~wasObserved, continue; end
        ind_j = ind_i & dates(:,1) == yr;
        dat = data(ind_j,:); % subset data for vessel i, year j
        nd = height(dat);
        if nd == 1, continue; end % A single AIS record is insufficient to determine visit duration
        gaps = diff(dat.move_id) ~= 1; % index when vessels leave the Southern Ocean
        gaps_ind = find(gaps);
        gaps_start = [1; gaps_ind + 1];
        gaps_end = [gaps_ind; nd];
        % Omit visits with a single AIS ping as these data cannot be used 
        % to determine visit duration
        singlePing = gaps_start == gaps_end;
        if any(singlePing)
            dat(gaps_start(singlePing),:) = [];
            nd = height(dat);
            gaps = diff(dat.move_id) ~= 1; % index when vessels leave the Southern Ocean
            gaps_ind = find(gaps);
            gaps_start = [1; gaps_ind + 1];
            gaps_end = [gaps_ind; nd];
        end
        alltimes = datenum(dat.date_time); % time of each AIS ping
        dlon = dat.longitude;
        dlat = dat.latitude;
        which_lon = lonmin < dlon & dlon <= lonmax;
        which_lat = latmin < dlat & dlat <= latmax;
        which_lon = repmat(which_lon, 1, 1, nlat);
        which_lat = repmat(reshape(which_lat, nd, 1, nlat), 1, nlon, 1);
        cell_index_ = repmat(reshape(cell_index, 1, nlon, nlat), nd, 1, 1);
        where = cell_index_(which_lon & which_lat); % the grid cell occupied during each AIS ping
        nvisits = sum(gaps) + 1;
        moves = diff(where) ~= 0 | gaps; % index when vessels move between grid cells
        moves_ind = find(moves);
        moves_start = [1; moves_ind + 1];
        moves_end = [moves_ind; nd];
        moves = [moves_start moves_end];
        v = 0;
        while v < nvisits
            v = v + 1;
            visit_ind = false(nd, 1);
            visit_ind(gaps_start(v):gaps_end(v)) = true;
            moves_v = moves(moves(:,1) == gaps_start(v) | ...
                moves(:,2) == gaps_end(v),:);
            nplaces = size(moves_v, 1);
            p = 0;
            while p < nplaces
                p = p + 1;
                moves_p = moves_v(p,:);
                place_ind = false(nd, 1);
                place_ind(moves_p(1):moves_p(2)) = true;
                alltimes_p = alltimes(place_ind);
                visitTime = alltimes_p(end) - alltimes_p(1); % time within grid cell
                % Account for times taken to move between grid cells
                if 1 < p && p < nplaces
                    visitTime = visitTime + ...
                        0.5 * (alltimes_p(1) - alltimes(moves_p(1) - 1)) + ...
                        0.5 * (alltimes(moves_p(2) + 1) - alltimes_p(end));
                elseif p == 1 && nplaces > 1
                    visitTime = visitTime + ...
                        0.5 * (alltimes(moves_p(2) + 1) - alltimes_p(end)); % time to move to other cell
                elseif p == nplaces && nplaces > 1
                    visitTime = visitTime + ...
                        0.5 * (alltimes_p(1) - alltimes(moves_p(1) - 1));
                end
                where_ = unique(where(place_ind));
                if length(where_) > 1, error('Error: location indexing problem. We need to index a single grid cell but multiple are selected.'); end
                shipTime_gridded(i,j,where_) = shipTime_gridded(i,j,where_) + visitTime;
            end
        end
    end
end


% Create a table of results to save as a csv file
vessel_id = repmat(all_vessel_id, 1, nyrs, ncells);
vessel_id = vessel_id(:);

vessel_name_u = unique(data(:,["vessel_id", "vessel_name"]), 'stable');
vessel_name_tmp = table(vessel_id);
vessel_name_tmp = join(vessel_name_tmp, vessel_name_u);
vessel_name = vessel_name_tmp.vessel_name;


year = repmat(reshape(yrs, 1, nyrs), nvessels, 1, ncells);
year = year(:);

lonmin_ = repmat(reshape(lonmin, 1, 1, nlon), nvessels, nyrs, 1, nlat);
lonmin = lonmin_(:);
lonmax_ = repmat(reshape(lonmax, 1, 1, nlon), nvessels, nyrs, 1, nlat);
lonmax = lonmax_(:);
latmin_ = repmat(reshape(latmin, 1, 1, 1, nlat), nvessels, nyrs, nlon, 1);
latmin = latmin_(:);
latmax_ = repmat(reshape(latmax, 1, 1, 1, nlat), nvessels, nyrs, nlon, 1);
latmax = latmax_(:);

ship_time = shipTime_gridded(:);

output_3x1 = table(vessel_id, vessel_name, year, lonmin, lonmax, latmin, latmax, ship_time);





switch saveOutput, case true
    path = fullfile(baseDirectory, 'MatLab', 'temp');
    if ~exist(path, 'dir')
        mkdir(fileparts(path), 'temp')
        addpath(genpath(path))
    end

    filename = 'ship_time_res_9x3.csv';
    filepath = fullfile(path, filename);
    writetable(output_9x3, filepath)

    filename = 'ship_time_res_3x1.csv';
    filepath = fullfile(path, filename);
    writetable(output_3x1, filepath)

end
















