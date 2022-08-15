% Download GLORYS model output directly from MatLab

% For instructions see:
% https://help.marine.copernicus.eu/en/articles/4799385-can-i-download-copernicus-marine-data-via-r-or-matlab#h_db0930094c

% This looks like a wrapper for the Python command line method for
% downloading the data. It should still be easier to use than Python, but
% will have all the same dependencies...

% Copernicus Marine Credentials 
user = 'ahunter';     %insert your username
pwd = 'CMEMS_Aidan_22';      %insert your password

%Download parameters
motu_server = 'https://my.cmems-du.eu/motu-web/Motu'; %server
productID = 'GLOBAL_MULTIYEAR_PHY_001_030';
datasetID = 'cmems_mod_glo_phy_my_0.083_P1M-m';

% Directories
thisFile = which('downloadGLORYS.m');
project = 'CUPIDO-risk-map';
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
baseDirectory = strrep(baseDirectory, ' ', '\ '); % account for spaces, Linux format

% addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))
% fullfile(baseDirectory, 'data')

%output directory to save the file
out_dir = fullfile(baseDirectory, 'data', 'physical_models', 'Copernicus_Programme', 'Mercator_Ocean_International', 'GLORYS', 'Southern\ Ocean');

% Choose area and variables
lon = ["-180", "179.9167"]; % circumpolar
lat = ["-80", "-45"];       % entire Southern Ocean
depth = ["0.494", "5727.9170"]; % all modelled depths
variables = ["thetao", "so"]; % temperature & salinity

% Download data one month at a time to satisfy download size limits
years = string(1993:2020); % get all available years
nyrs = length(years);
months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"];
days = ["15", "16"]; % day is 15 for Feb and 16 for all other months
times = ["12:00:00", "00:00:00"];
niter = (nyrs-1) * 12 + 5; % data available for all months from 1993-2019, and for 1st 5 months of 2020

pythonVersion = '3';

showProgressBar = true;

% Download data in chunks -- this will take a long time and require lots of memory
for i = 1:nyrs
    year = years(i);
    for j = 1:length(months)
        fileCount = (i - 1) * 12 + j; % keep track of how many data files have been processed
        if fileCount == 1
            switch showProgressBar, case true
                progress = waitbar(0, 'Calculating...');
            end
        end
        month = months(j);
        if strcmp(years(i), "2020") && str2double(month) > 5, break; end
        isFeb = strcmp(month, '02');
        if isFeb, day = days(1); else, day = days(2); end
        filename = append(year, month, day, ".nc");
        fileExists = exist(filename, 'file') == 2;
        switch fileExists, case true
            switch showProgressBar, case true
                waitbar(fileCount / niter, progress)
            end
            continue
        end
        for k = 1:length(times)
            time = times(k);
            date = repmat(append(year, "-", month, "-", day, " ", time), 1, 2);
            %%Build the command line for motuclient
            motu_line = append("python", pythonVersion, " -m motuclient --motu ", motu_server," --service-id ", productID, "-TDS --product-id ", datasetID, " --longitude-min ", lon(1), " --longitude-max ", lon(2), " --latitude-min ", lat(1), " --latitude-max ", lat(2), " --date-min ", char(34), date(1), char(34), " --date-max ", char(34), date(2), char(34), " --depth-min ", depth(1), " --depth-max ", depth(2), " --variable ", variables(1), " --variable ", variables(2), " --out-dir ", out_dir, " --out-name ", filename, " --user ", user, " --pwd ", pwd);
            %%Run the command line
            system(motu_line)
        end
        switch showProgressBar, case true
            waitbar(fileCount / niter, progress)
        end
        pause(0.25)
    end
end



