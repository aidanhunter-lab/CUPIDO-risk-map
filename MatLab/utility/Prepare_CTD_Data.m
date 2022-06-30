function [Data, MetaData, cruises] = Prepare_CTD_Data(project, source, region, varargin)
% Load and organise CTD data sampled from Region and downloaded from Source.
% This is a wrapper function.

%% Preamble -- extract/define optional arguments

extractVarargin(varargin)

% Extra arguments for loading data
if ~exist('displayData', 'var'), displayData = true; end
if ~exist('printCruiseStructure', 'var'), printCruiseStructure = false; end % true => structure containing cruise metadata is printed on screen
if ~exist('printFileStructure', 'var'), printFileStructure = false; end % true => structure containing data file names is displayed on screen

% Extra arguments for cleaning data
if ~exist('plotMaps', 'var'), plotMaps = true; end
if ~exist('flagGoodMeasure', 'var'), flagGoodMeasure = [49 50 56]; end % examine NetCDF files to check measurement quality control flags
if ~exist('plotSpuriousPoints', 'var'), plotSpuriousPoints = true; end
if ~exist('removeSpuriousPoints', 'var'), removeSpuriousPoints = true; end
if ~exist('onlyDepthProfiles', 'var'), onlyDepthProfiles = false; end
% parameter used to identify depth profiles
if ~exist('CoV_lim', 'var'), CoV_lim = 0.05; end
% parameters used to (loosely) identify spurious measurements
if ~exist('cutOffValue', 'var'), cutOffValue = 0.5; end                    % test moving std against cutOffValue
if ~exist('leeway', 'var'), leeway = 2; end                              % non-negative integer -- extra data points to test before specifying a point as spurious.
if ~exist('movstdRange', 'var'), movstdRange = 0.1; end                    % range of moving-std as fraction of total points
if ~exist('XaxisLocation', 'var'), XaxisLocation = 'bottom'; end
if ~exist('ShowMapTitle', 'var'), ShowMapTitle = true; end


%% Load the data

[Data, MetaData, cruises] = Load_CTD_Data(project, source, region, ...
    'displayData', displayData, 'printCrusieStructure', printCruiseStructure, ...
    'printFileStructure', printFileStructure);


%% Clean the data

[Data, MetaData, cruises] = Clean_CTD_Data(project, source, region, ...
    Data, MetaData, cruises, ...
    'plotMaps', plotMaps, 'flagGoodMeasure', flagGoodMeasure, 'plotSpuriousPoints', plotSpuriousPoints, ...
    'removeSpuriousPoints', removeSpuriousPoints, 'onlyDepthProfiles',onlyDepthProfiles, ...
    'CoV_lim', CoV_lim, 'cutOffValue', cutOffValue, 'leeway', leeway, 'movstdRange', movstdRange, ...
    'XaxisLocation', XaxisLocation, 'ShowMapTitle', ShowMapTitle);

