%% Generate risk map

% Build up a risk map gradually layer by layer. Focus on a single year...


%% Directories
% Adjust search path to include all MatLab scripts and the 'data' directory
project = 'CUPIDO-risk-map';
thisFile = which('run_risk_map.m');
baseDirectory = thisFile(1:strfind(thisFile, project)+length(project)-1);
addpath(genpath(fullfile(baseDirectory, 'MatLab')))
addpath(genpath(fullfile(baseDirectory, 'data')))

%% Model domain
% Specify the modelled area and horizontal & vertical resolutions
domain = map_domain();
% Returns info on domain grid resolution [lon-lat coords] and grid cell
% areas & volumes [m^2 & m^3].


%% Model parameters
% The sinking speed equations are written in cgs units.
pars = initialise_parameters(domain);
% - fpCm [g C / cm^3], faecal pellet mean carbon content
% - FPprod [mm^3 / h / g dry mass], faecal pellet production rate per krill dry mass

%% Prepare model forcing data
[forc, domain] = prepare_forcing(baseDirectory, domain, 'loadFromFile', true);
% Returns model forcing data:
% - Krill [individuals/m^3]
% - Krill mean length [mm] on log_e scale
% - Krill SD length [mm] on log_e scale
% - Krill mean dry weight [mg] on log_e scale
% - Krill SD dry weight [mg] on log_e scale
% - Krill mean wet weight [mg] on log_e scale
% - Krill SD wet weight [mg] on log_e scale
% - Krill mean C [mg] on log_e scale
% - Krill SD C [mg] on log_e scale

% THE LENGTH DATA MAY NEED TO BE HANDLED DIFFERENTLY.
% MAYBE BEST TO NOT RESOLVE LENGTH DATA BY MAP GRID CELL BECAUSE AVERAGE
% KRILL SIZE IN THE SWARM WILL BE SOMEWHAT RANDOM. SUCH A LARGE RANGE IN
% MEAN LENGTHS ACROSS GRID CELLS WILL RESULT IN VERY DIFFERENT FAECAL
% PELLET PRODUCTION RATES. PERHAPS THIS IS NOT SUCH A BAD THING SINCE IT
% PROVIDES A NATURAL WAY TO EXAMINE THE EFFECT OF VARYING KRILL SIZE; WE
% CAN JUST SEPARATELY EXAMINE GRID CELLS WITH DIFFERENT SIZED KRILL.


%% Initial values
% Generate initial values for all modelled variables.
init = generate_initials(domain, pars, forc);
% - CFP [g C in faecal pellets / grid cell]


%% Integrate equations

% Set options for differential equation solver
[odeIntegrator, odeOptions] = integration_options(pars);

% Parallelise integrations over mapped grid cells
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Integrate to solve equations
tic; fprintf('\n'); disp(append("started: ", string(datetime('now'))))
[out, outStruct] = integrateModel(domain, pars, forc, init, odeIntegrator, odeOptions);
integrationTime = toc / 60; fprintf('\n'); disp(append("finished: ", string(datetime('now')))); fprintf('\n');
disp(['Integration time: ' num2str(floor(integrationTime)) ' mins, ' ...
    num2str(floor(mod(60*integrationTime,60))) ' secs'])

% PERHAPS THE SINKING SPEED EQUATIONS COULD BE MODIFIED TO ACCOUNT FOR
% VOLUME REDUCTIONS RESULTING FROM REMINERALISATION. I DON'T WANT TO MODEL
% FAECAL PELLETS OF VARYING SIZE (THEY WILL BE CONSTANT VOLUME), HOWEVER, I
% DO WANT TO MODEL REMINERALISATION WHICH IS A PROCESS THAT REDUCES SIZE
% AND THEREFORE REDUCES SINK SPEED. REMINERALISATION WILL BE MODELLED AS
% REDUCTION IN TOTAL CARBON (REDUCING THE NUMBER OF FAECAL PELLETS, BUT NOT
% THEIR SIZE). I COULD TRY TO MODIFY THE SINK SPEED EQUATION TO ACCOUNT FOR
% THIS BY REDUCING SINK SPEED WITHIN DEEPER LAYERS WHERE PELLETS HAVE SPENT
% MORE TIME BEING DEGRADED. THIS WILL TAKE A BIT OF THOUGHT, BUT IS WORTH
% KEEPING IN MIND... I'M NOT SURE HOW LARGE/IMPORTANT THE EFFECT WILL BE...


%% Include plastics

% Plastic dimensions. These may be comparable to the faecal pellet
% dimensions. Look up a paper that shows plastic particles inside faecal
% pellet and get approx dimensions from there.
% The plastics that are excreted may be smaller than those that are eaten
% because they are ground down to smaller sizes. An ingested piece of
% plastic may not be excreted all at once, it be excreted gradually as
% smaller pieces.

% Vp = 229972525; % mu m^3
% 1e-12 * Vp




%% Calculate faecal pellet production rates


















