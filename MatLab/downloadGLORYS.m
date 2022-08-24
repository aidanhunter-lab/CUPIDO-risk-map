%% Download GLORYS model output directly from MatLab

% For instructions and dependencies see:
% https://help.marine.copernicus.eu/en/articles/4799385-can-i-download-copernicus-marine-data-via-r-or-matlab#h_db0930094c

%% Copernicus Marine Credentials

% Insert your username
user = 'ahunter';
% and your password
pwd = 'CMEMS_Aidan_22';

%% Download parameters for GLORYS model

% Server
motu_server = 'https://my.cmems-du.eu/motu-web/Motu';
% Model type
productID = 'GLOBAL_MULTIYEAR_PHY_001_030';
% GLORYS data are available with temporal resolution as a daily or monthly
% average
timeResolutions = {'monthly', 'daily'};
% Choose time resolution for download. Select 'monthly' for data informing
% the risk map, and select 'daily' for validating GLORYS model against
% available CTD data.
timeRes = timeResolutions{2};
% Model details
switch timeRes
    case 'monthly'
        datasetID = 'cmems_mod_glo_phy_my_0.083_P1M-m';
    case 'daily'
        datasetID = 'cmems_mod_glo_phy_my_0.083_P1D-m';
end

%% Directories

% thisFile = which('downloadGLORYS.m');
% project = 'CUPIDO-risk-map';
% baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
% baseDirectory = strrep(baseDirectory, ' ', '\ '); % account for spaces, Linux format

% AS FILE SIZES ARE LARGE AND WE'LL DOWNLOAD LOTS OF THEM, SELECT A STORAGE
% DIRECTORY OUTSIDE OF THE GIT REPO. LATER WE WILL COPY REQUIRED FILES INTO
% THE GIT REPO, AS NEEDED, AND STORE THEM WITH GIT-LFS. IF TOO MANY OF
% THESE LARGE FILES ARE STORED IN THE GIT REPO THEN THE REPO SLOWS DOWN AND
% BECOMES FRUSTRATING TO USE.

% Directory to store GLORYS data -- outside of Git Repo
dataDirectory = ['/home/aihunt/Documents/work/CUPIDO/data/physical_models/' ...
    'Copernicus_Programme/Mercator_Ocean_International/GLORYS'];
switch timeRes
    case 'monthly'
        subDirectory = 'Southern Ocean';
    case 'daily'
        subDirectory = 'Scotia Sea/2003 Jan-Feb';
end
out_dir = fullfile(dataDirectory, subDirectory);
% If directory does not exist then create it
dir_exists = exist(out_dir, 'dir') == 7;
if ~dir_exists, mkdir(dataDirectory, subDirectory); end

% Code spaces with a backslash (Linux style)
out_dir = strrep(out_dir, ' ', '\ ');


%% Choose area and variables

% Select depth range -- all available depths
depth = ["0.494", "5727.9170"];
% Download temperature & salinity variables
variables = ["thetao", "so"];

% Select grid cells to download -- this depends on the chosen time
% resolution. We want the entire Southern Ocean (monthly data) for the risk
% map, and the Scotia Sea (daily data) for model validation.
switch timeRes
    case 'monthly'
        lon = ["-180", "179.9167"]; % circumpolar
        lat = ["-80", "-45"];       % entire Southern Ocean
    case 'daily'
        lon = ["-63", "-26"];
        lat = ["-63", "-53"];
end


%% Download the data

pythonVersion = '3';

switch timeRes
    case 'monthly'
        % Download data one month at a time to satisfy download size limits
        years = string(1993:2020); % get all available years
        nyrs = length(years);
        months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"];
        days = ["15", "16"]; % day is 15 for Feb and 16 for all other months
        times = ["12:00:00", "00:00:00"];
        niter = (nyrs-1) * 12 + 5; % data available for all months from 1993-2019, and for 1ast 5 months of 2020

        showProgressBar = true;

        % Download data in chunks -- this will take a long time and require
        % lots of memory (approx 0.5 GB per file). If too much memory is
        % required then download data from only the most recent years.
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
                filepath = fullfile(out_dir, filename);
                fileExists = exist("filepath", 'file') == 2;
%                 fileExists = exist(filename, 'file') == 2;
                switch fileExists, case true
                    switch showProgressBar, case true
                        waitbar(fileCount / niter, progress)
                    end
                    continue
                end
                for k = 1:length(times)
                    time = times(k);
                    date = repmat(append(year, "-", month, "-", day, " ", time), 1, 2);
                    % Build the command line for motuclient
                    motu_line = append("python", pythonVersion, " -m motuclient --motu ", motu_server," --service-id ", productID, "-TDS --product-id ", datasetID, " --longitude-min ", lon(1), " --longitude-max ", lon(2), " --latitude-min ", lat(1), " --latitude-max ", lat(2), " --date-min ", char(34), date(1), char(34), " --date-max ", char(34), date(2), char(34), " --depth-min ", depth(1), " --depth-max ", depth(2), " --variable ", variables(1), " --variable ", variables(2), " --out-dir ", out_dir, " --out-name ", filename, " --user ", user, " --pwd ", pwd);
                    % Run the command line
                    system(motu_line)
                end
                switch showProgressBar, case true
                    waitbar(fileCount / niter, progress)
                end
                pause(0.25)
            end
        end

    case 'daily'

        year = "2003"; % I have CTD data from 2003
        months = ["01", "02"]; % Jan and Feb
        days = {num2str((7:31)', '%02i'), num2str((1:17)', '%02i')};
        time = "12:00:00";
        date = {append(year, "-", months(1), "-", days{1}(1,:), " ", time), ...
            append(year, "-", months(2), "-", days{2}(end,:), " ", time)};

        filename = append(append(year, months(1), days{1}(1,:)), "-", ...
            append(year, months(2), days{2}(end,:)), "_daily", ".nc");
        filepath = fullfile(out_dir, filename);
        fileExists = exist("filepath", 'file') == 2;

        % The daily data that I am using can be downloaded as a single
        % file, but if you need a larger file that exceeds the download
        % limit then write code to download in chunks -- similar to the
        % monthly case above.

        % Build the command line for motuclient
        motu_line = append("python", pythonVersion, " -m motuclient --motu ", motu_server," --service-id ", productID, "-TDS --product-id ", datasetID, " --longitude-min ", lon(1), " --longitude-max ", lon(2), " --latitude-min ", lat(1), " --latitude-max ", lat(2), " --date-min ", date{1}, " --date-max ", date{2}, " --depth-min ", depth(1), " --depth-max ", depth(2), " --variable ", variables(1), " --variable ", variables(2), " --out-dir ", out_dir, " --out-name ", filename, " --user ", user, " --pwd ", pwd);
        % Run the command line
        system(motu_line)

end
