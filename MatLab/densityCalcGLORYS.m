%% Calculate water density for GLORYS model output.

% This script will take a long time to complete, but we should only need to
% run it once.

% The monthly average GLORYS output should already be downloaded by running
% downloadGLORYS.m. 

% All of the monthly data for the Southern Ocean should be stored in a
% sub-directory of the data directory. These data are separate files for
% each month in each year (1993-2020) and the files are named in the 
% year-month-day format, e.g., 20100516.nc corresponds to the 16th of May 
% 2010. The data are monthly averages all recorded as either day 15 or 16
% of the month.

% There are 329 separate data sets, each approximately 360MB. All 12 months
% of 1993-2019, and months 1-5 of 2020. After running this script, each
% file will increase in size to 540MB.


%% Preamble

% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('densityCalcGLORYS.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))
dataPath = fullfile(baseDirectory, 'data', 'physical_models', ... 
    'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', ...
    'Southern Ocean');
months = ["01","02","03","04","05","06","07","08","09","10","11","12"];
days = ["15","16"];
% Choose a year range
ymin = 1993;
ymax = 2020;
yrs = ymin:ymax;
nyrs = length(yrs);
ntimes = 12 * (nyrs - 1) + 5; % number of data sets (see above description)


%% Density calculation

% Iterate the density calculations, one data set at a time. This is time
% consuming, but necessary as the data are large.
% With each loop, a data is loaded, the density variable is calculated,
% then the new variable is saved into the data set for future use.

showProgressBar = true;

for i = 1:nyrs
    nmonths = 12;
    if i == nyrs, nmonths = 5; end
    for j = 1:nmonths
        fileCount = (i - 1) * 12 + j; % keep track of how many data files have been processed
        f = {fullfile(dataPath, append(num2str(yrs(i)), months(j), days(1), ".nc")), ...
            fullfile(dataPath, append(num2str(yrs(i)), months(j), days(2), ".nc"))}; % possible file names for this year & month
        filename = f{[exist(f{1}, 'file'), exist(f{2}, 'file')] == 2}; % choose the file name that exists in the search path
        if fileCount == 1
            switch showProgressBar, case true
                progress = waitbar(0, 'Calculating...');
            end
        end
%             ncdisp(filename)
        fileInfo = ncinfo(filename);
        % If the density variable has already been calculated then skip to
        % next loop iteration -- useful if this process is done
        % incrementally, which is likely because it takes a while.
        skip = any(strcmp({fileInfo.Variables.Name}, 'density'));
        if skip
            switch showProgressBar, case true
                waitbar(fileCount / ntimes, progress)
            end
            continue
        end
        % Get data dimensions
        dataDim = [{fileInfo.Dimensions.Name};
            {fileInfo.Dimensions.Length}];
        nlon = dataDim{2,strcmp(dataDim(1,:), 'longitude')};
        nlat = dataDim{2,strcmp(dataDim(1,:), 'latitude')};
        ndep = dataDim{2,strcmp(dataDim(1,:), 'depth')};
        % Load NetCDF data
        dat.depth = ncread(filename, 'depth');
        dat.longitude = ncread(filename, 'longitude');
        dat.latitude = ncread(filename, 'latitude');
        dat.temp = ncread(filename, 'thetao');
        dat.salin = ncread(filename, 'so');
        % Calculate water pressure from depth & latitude
        dat.pressure = gsw_p_from_z(repmat(reshape(-dat.depth, ...
            [1 ndep]), [nlat 1]), dat.latitude);
        % Convert from practical salinity to absolute salinity using gsw_SA_from_SP.m
        dat.temp = permute(dat.temp, [3 1 2]);
        dat.salin = permute(dat.salin, [3 1 2]);
        lon3D = repmat(reshape(dat.longitude, [1 nlon]), [ndep 1 nlat]);
        lat3D = repmat(reshape(dat.latitude, [1 1 nlat]), [ndep nlon 1]);
        press3D = repmat(reshape(permute(dat.pressure, [2 1]), ...
            [ndep 1 nlat]), [1 nlon 1]);
        dat.temp = permute(dat.temp(:,:), [2 1]);
        dat.salin = permute(dat.salin(:,:), [2 1]);
        press3D = permute(press3D(:,:), [2 1]);
        lon3D = permute(lon3D(:,:), [2 1]);
        lat3D = permute(lat3D(:,:), [2 1]);
        % gsw_SA_from_SP.m is a memory-hungry function so call it
        % iteratively for each depth layer
        for k = 1:ndep
            dat.salin(:,k) = gsw_SA_from_SP(dat.salin(:,k), press3D(:,k), ...
                lon3D(:,k), lat3D(:,k));
        end
        % Convert potential temperature into conservative temperature
        dat.temp = gsw_CT_from_pt(dat.salin, dat.temp);
        % Calculate water density from temperature, salinity and pressure
        dat.density = gsw_rho(dat.salin, dat.temp, press3D);
        dat.density = permute(reshape(permute(dat.density, [2 1]), [ndep nlon nlat]), [2 3 1]);
        % For storage efficiency, minimise memory usage by packaging the
        % density variable as an integer, unpacked using scale and offset
        % parameters.
        % Useful instructions found here: https://www.unidata.ucar.edu/software/netcdf/workshops/2010/bestpractices/Packing.html
        DataType = 'int16';
        add_offset = double(min(dat.density(:)));
        prec = 7e2; % Chosen to be compatible with the int16 storage
        scale_factor = double(1 / prec);
        FillValue = -32767; % NaN values represented by scalar when int16 is used to store data
        dat.density_packed = (dat.density - add_offset) / scale_factor;
        dat.density_packed(isnan(dat.density)) = FillValue;
        dat.density_packed = int16(dat.density_packed);
        % Create new density variable in the NetCDF file
        nccreate(filename, 'density', 'Dimensions', ... 
            {'longitude', nlon, 'latitude', nlat, 'depth', ndep, 'time', 1}, ...
            'Datatype', DataType) % , 'FillValue', FillValue)
        % Save the (integer-transformed) density variable into the NetCDF file
        ncwrite(filename, 'density', dat.density_packed)
        % Store data attributes -- in same format as temperature and salinity
        ncwriteatt(filename, 'density', '_FillValue', FillValue)
        ncwriteatt(filename, 'density', 'long_name', 'Density')
        ncwriteatt(filename, 'density', 'standard_name', 'sea_water_density')
        ncwriteatt(filename, 'density', 'units', 'kg/m^3')
        ncwriteatt(filename, 'density', 'unit_long', 'kilograms per metre cubed')
        ncwriteatt(filename, 'density', 'add_offset', add_offset)
        ncwriteatt(filename, 'density', 'scale_factor', scale_factor)
        clearvars dat lon3D lat3D press3D
        switch showProgressBar, case true
            waitbar(fileCount / ntimes, progress)
            pause(0.25)
        end
    end
end

