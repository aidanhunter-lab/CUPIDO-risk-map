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


%% Prepare model forcing data
[forc, domain] = prepare_forcing(baseDirectory, domain);


%% Model parameters
% The sinking speed equations are written in cgs units.
pars = initialise_parameters(domain);


%% Initial values
% Generate initial values for all modelled variables.
% The modelled variable is total carbon in krill faecal pellets.
init = generate_initials(domain, pars, forc);


%% Integrate equations

% Set solver options.
[odeIntegrator, odeOptions] = integration_options(pars);

out = integrateModel(domain, pars, forc, init, odeIntegrator, odeOptions);






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

Vp = 229972525; % mu m^3
1e-12 * Vp




%% Calculate faecal pellet production rates


















