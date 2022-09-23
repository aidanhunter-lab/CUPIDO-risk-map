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
% The integrator is called separately for each distinct period of forcing
% data (which is monthly intervals).
odeMaxTime = pars.dt_max; % max integration timestep (seconds)
odeInitTime = 0.5 * odeMaxTime; % initial integration timestep (solver will automatically reduce this if required)
odeOptions = odeset('InitialStep', odeInitTime, 'MaxStep', odeMaxTime); % Integration tolerances can be set here if required...
% Solver functions
integratorChoices = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s'};
odeIntegrator = integratorChoices{1};
odeSolve = str2func(odeIntegrator);

% WRITE A FUNCTION TO INTEGRATE THE MODEL. THIS WILL LOOP THROUGH THE
% FORCING DATA MONTHS, CALLING THE INTEGRATOR SEPARATELY FOR EACH PERIOD. I
% CAN ALSO USE PARALLEL PROCESSING TO QUICKLY HANDLE THE HORIZONTAL GRID
% CELLS AS THE MODEL IS IDENTICALLY RUN FOR EACH UNIQUE CELL.
out = integrateModel(domain, pars, forc, init, odeIntegrator, odeOptions);



sol = odeSolve(@(t, y) risk_map_model(t, y, domain, pars, forc), [0 1], init.CFP, odeOptions);

% sol = odeSolve(@(t, v_in) ODEs(t, v_in, parameterList, forcing, j, false), [0 1], v_in, odeOptions);


% Integrate
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ...
    odeIntegrator, odeOptions, 'returnExtra', returnExtra);


%% Calculate sinking speeds

% Find faecal pellet sinking speed, cm / s
[v, Re] = sinking_speed(pars.shape, 1e-3 * dat.density, pars.rho_p, pars.mu, pars.g, ...
    'L', pars.L, 'D', pars.D, 'Dl', pars.Dl, 'Di', pars.Di, 'Ds', pars.Ds, ...
    'useMeans_Re', false, 'returnMax_Re', true);

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


















