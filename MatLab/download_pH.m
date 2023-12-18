%% Download pH data directly from MatLab

% This requires Python is installed
pythonVersion = '3';

% For instructions and dependencies see:
% https://help.marine.copernicus.eu/en/articles/4799385-can-i-download-copernicus-marine-data-via-r-or-matlab#h_db0930094c

%% Copernicus Marine Credentials

% Store your Copernicus username & password in a hidden text file called
% '.Copernicus Marine Credentials.txt', formatted as a table with column
% names 'user' & 'pwd'.
cmc = readtable('.Copernicus Marine Credentials.txt');

% Insert your username
user = cmc.user {1};
% and your password
pwd = cmc.pwd{1};

%% Download parameters for ESA SST CCI and C3S reprocessed sea surface temperature analyses

% Download server
motu_server = 'https://my.cmems-du.eu/motu-web/Motu';
% Model services
serviceID = 'MULTIOBS_GLO_BIO_CARBON_SURFACE_REP_015_008-TDS';
% Model product
productID = 'dataset-carbon-rep-monthly';


%% Directories
thisFile = which('download_pH.m');
dirBase = fileparts(fileparts(thisFile));
dirData = fullfile(dirBase, 'data', 'pH', 'SOCAT');
addpath(genpath(dirData))

out_dir = fullfile(dirData);
% If directory does not exist then create it
dir_exists = exist(out_dir, 'dir') == 7;
if ~dir_exists, mkdir(dirData); end

% Code spaces with a backslash (Linux style)
out_dir = strrep(out_dir, ' ', '\ ');


%% Choose download variables

variables = ["omega_ar", "omega_ar_uncertainty", "omega_ca", "omega_ca_uncertainty", ...
    "ph", "ph_uncertainty"];

%% Download the data

% These are monthly data with a 1x1 degree (lat-lon) resolution.

% Build the command line for motuclient.
% Output variables
vars_ = arrayfun(@(z) sprintf('%s ', '--variable', z), variables, 'UniformOutput', false);
vars = ' ';
for i = 1:length(variables), vars = append(vars, vars_{i}); end
vars = vars(1:end-1);

% Coordinates
lonrange = [-180, 180];
latrange = [-90, -45];
lonrange_ = string(lonrange);
latrange_ = string(latrange);
reslon = 1;
reslat = 1;
longrid = lonrange(1):reslon:lonrange(2);
latgrid = latrange(1):reslat:latrange(2);
nlon = -1 + length(longrid);
nlat = -1 + length(latgrid);

yrrange = [1985, 2021];

filename = append("pH_", string(yrrange(1)), "0115_", string(yrrange(2)), "1215.nc");

% Generate & store the Python code to download the daily data
code = append("python", pythonVersion, " -m motuclient --motu ", ...
        motu_server, " --service-id ", serviceID, " --product-id ", productID, ...
        " --longitude-min ", lonrange_(1), " --longitude-max ", lonrange_(2), ...
        " --latitude-min ", latrange_(1), " --latitude-max ", latrange_(2), ...
        " --date-min ", string(yrrange(1)), "-01-15 00:00:00", " --date-max ", ...
        string(yrrange(2)), "-12-15 00:00:00", vars, " --out-dir ", out_dir, ...
        " --out-name ", filename, " --user ", user, " --pwd ", pwd);

% Download
system(code)

