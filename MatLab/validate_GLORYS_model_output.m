% Validate GLORYS physical model data using CTD samples

%% Preamble

% Adjust search path to include all MatLab scripts and the 'data' directory
addpath(genpath(fileparts(which('validate_GLORYS_model_output.m'))))
addpath(genpath(fullfile(fileparts(fileparts(which('validate_GLORYS_model_output.m'))), 'data')))


%% Load & clean CTD data

% Set project, data source, and region -- in accordance with directory structure
project = 'CUPIDO-risk-map';
source = 'BODC';
region = 'South Georgia';

[Data, MetaData, ~] = Prepare_CTD_Data(project, source, region, ... 
    'plotMaps', true, 'displayData', false, ...
    'flagGoodMeasure', [49 50 56], 'plotSpuriousPoints', false, ... 
    'removeSpuriousPoints', true, 'onlyDepthProfiles', true, ...
    'ShowMapTitle', true);


% To validate the GLORYS model output let's use CTD samples from cruises
% JR20141115 (2014-15), JR15002 (2015-16), and JR16003 (2016-17).

cruises = {'JR20141115', 'JR15002', 'JR16003'};

% For each cruise we need to find the range of times and locations sampled,
% then match these to the GLORYS model output. To avoid downloading/storing
% excessivley large files, we only retain model output matching the CTD
% samples.

% The BODC data has time variable corresponding to days since 1/1/-4713,
% which differs from MatLab's default of 1/1/0000.
BODC_startDate = [-4712 1 1 0 0 0];
adjustTime = datenum(BODC_startDate);

% For each cruise, note the GLORYS output data to extract...
for i = 1:length(cruises)
    % Subset data by cruise
    d = MetaData(strcmp(MetaData.Cruise, cruises(i)),:);
    timeRange = [min(d.Time) max(d.Time)] + adjustTime;
%     dateRange = datevec(timeRange);
    dateRange2 = datestr(timeRange);
    lonRange = [min(d.Longitude) max(d.Longitude)];
    latRange = [min(d.Latitude) max(d.Latitude)];
    extractData.(cruises{i}).dates = dateRange2;
    extractData.(cruises{i}).lons = lonRange;
    extractData.(cruises{i}).lats = latRange;
end














