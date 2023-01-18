%% Estimate microplastic emissions from Antarctic facilities
% Use info on population sizes at facilities, combined with assumptions
% about MP emission rates per person.

%% Preamble
% Include all required (sub)directories within the search path
project = 'CUPIDO-risk-map';
thisFile = which('plastic_emissions_Antarctic_facilities.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))

%% Load/organise data

filename = 'COMNAP_Antarctic_Facilities_Master.csv';
filepath = fullfile(baseDirectory, 'data/research stations/COMNAP');
stations = readtable(fullfile(filepath, filename));

% Omit rows lacking info on population size
stations = stations(~isnan(stations.Peak_Population),:);

%% Calculate

% Keep everything as simple as possible to begin with...
% Assume that each facility has a population size equal to the peak
% population listed in the data, and that this population size is constant
% throughout time (or at least a single year). We can then find the
% person-days at each facility. Next, we estimate MP emissions per person
% by making assumptions about laundry, hygiene products. This will yield a
% basic estimate of MP emssions per facility. These estimates will be
% biased by (1) underestimating per-person emissions and (2) overestimated
% facility population size.

% Population size, person-days
pop = stations.Peak_Population; % assume this represents facility population size throughout a single year
pdays = 365 .* pop; % person-days (person day / year) at each facility

% Microbeads, hygiene products (Gouin, 2015)
mpbead_rate = 17.5; % MP bead typical usage rate [mg / day / person]
mpbead_rate_range = [7.5, 27.5];

% Laundry fibres (Napper 2016)
washesPerWeek = 1; % weekly rate at which a single person does laundry
washesPerDay = washesPerWeek / 7;
washSize = 6; % weight [kg] of material being laundered per wash
fibreDiameter_min = 1e-6 * 11.9; % m
fibreDiameter_max = 1e-6 * 17.7; % m
fibreLength_min = 1e-3 * 5.0; % m
fibreLength_max = 1e-3 * 7.8; % m
fibreVolume_min = fibreLength_min .* pi .* 1/4 .* fibreDiameter_min .^ 2; % m^3
fibreVolume_max = fibreLength_max .* pi .* 1/4 .* fibreDiameter_max .^ 2; % m^3

diameter_PolyesterCotton = 1e-6 * 17.74;
diameter_Polyester = 1e-6 * 11.91;
diameter_Acrylic = 1e-6 * 14.05;
length_PolyesterCotton = 1e-3 * 4.99;
length_Polyester = 1e-3 * 7.79;
length_Acrylic = 1e-3 * 5.44;
volume_PolyesterCotton = length_PolyesterCotton .* pi .* 1/4 .* diameter_PolyesterCotton .^ 2;
volume_Polyester = length_Polyester .* pi .* 1/4 .* diameter_Polyester .^ 2;
volume_Acrylic = length_Acrylic .* pi .* 1/4 .* diameter_Acrylic .^ 2;


densityNumeric_PolyesterCotton = 1e3 * 334800; % fibres / g
densityNumeric_Polyester = 1e3 * 475998; % fibres / g
densityNumeric_Acrylic = 1e3 * 763130; % fibres / g

density_PolyesterCotton = 1 / (volume_PolyesterCotton .* densityNumeric_PolyesterCotton); % g / m^3
density_Polyester = 1 / (volume_Polyester .* densityNumeric_Polyester); % g / m^3
density_Acrylic = 1 / (volume_Acrylic .* densityNumeric_Acrylic); % g / m^3

fibreReleaseRate_PolyesterCotton = 22992; % fibres / kg of polyester-cotton blend
fibreReleaseRate_Polyester = 82672; % fibres / kg of polyester
fibreReleaseRate_Acrylic = 121464; % fibres / kg of acrylic

% (washes / day / person) * (kg / wash) * (fibres / kg) = (fibres / day / person)
fibresPerDay_PolyesterCotton = washesPerDay .* washSize .* fibreReleaseRate_PolyesterCotton; % (fibres / day / person)
fibresPerDay_Polyester = washesPerDay .* washSize .* fibreReleaseRate_Polyester; % (fibres / day / person)
fibresPerDay_Acrylic = washesPerDay .* washSize .* fibreReleaseRate_Acrylic; % (fibres / day / person)

% Emission rate
MPrate_bead = mpbead_rate .* pdays; % mg MP / year

MPrate_fibres_PolyesterCotton = fibresPerDay_PolyesterCotton .* pdays; % (fibres / year)
MPrate_fibres_Polyester = fibresPerDay_Polyester .* pdays; % (fibres / year)
MPrate_fibres_Acrylic = fibresPerDay_Acrylic .* pdays; % (fibres / year)

MPrate_mass_fibres_PolyesterCotton = MPrate_fibres_PolyesterCotton ./ densityNumeric_PolyesterCotton; % (g / year)
MPrate_mass_fibres_Polyester = MPrate_fibres_Polyester ./ densityNumeric_Polyester; % (g / year)
MPrate_mass_fibres_Acrylic = MPrate_fibres_Acrylic ./ densityNumeric_Acrylic; % (g / year)


