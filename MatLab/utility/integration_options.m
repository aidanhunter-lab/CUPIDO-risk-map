function [odeIntegrator, odeOptions] = integration_options(pars, varargin)
% Specify integration function and parameters

extractVarargin(varargin)
integratorChoices = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s'};
if ~exist('integratingFunction', 'var')
    integratingFunction = 'ode45';
end
if ~ismember(integratingFunction, integratorChoices)
    error(['If optional argument "integratingFunction" is specified it must be one of' ...
        sprintf(' "%s",', integratorChoices{1:end-1}) ' or ' sprintf('"%s"', integratorChoices{end})])
end
odeIntegrator = integratingFunction;


if ~exist('odeMaxTime', 'var')
    odeMaxTime = pars.dt_max; % max integration timestep (seconds)
end
if odeMaxTime ~= pars.dt_max
    warning('Optional parameter "odeMaxTime" has been set to a value different from that specified in parameter struct "pars". Care must be taken here as model stability depends upon odeMaxTime!')
end
if ~exist('odeInitTime', 'var')
    odeInitTime = 0.5 * odeMaxTime; % initial integration timestep (the solver will automatically reduce/increase this if needed)
end
if odeInitTime > odeMaxTime
    odeInitTime = odeMaxTime;
    warning('Poor choice of optional parameters! "odeInitTime" has been reset to equal "odeMaxTime".')
end
if ~exist('RelTol', 'var')
    RelTol = 1e-3;
end
if ~exist('AbsTol', 'var')
    AbsTol = 1e-6;
end

odeOptions = odeset('InitialStep', odeInitTime, 'MaxStep', odeMaxTime, ...
    'RelTol', RelTol, 'AbsTol', AbsTol);




