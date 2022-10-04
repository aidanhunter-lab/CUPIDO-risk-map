function [out, outStruct] = integrateModel(domain, pars, forc, init, odeIntegrator, odeOptions)
% Calls an integrating function to solve equations within each modelled
% water column. As each water column in the map grid is characterised by 
% identical dynamics, use parallel processing to simultaneously solve
% equations across distinct water columns. To overcome discontinuities due
% to step-changes in forcing data, the integrating function must be called
% separately for each period of distinct forcing. We therefore integrate
% within a nested loop indexing over water column and forcing data period.


%% Preamble

nlon = domain.nlon;
nlat = domain.nlat;
ndepth = domain.ndepth;
nColumns = nlon * nlat; % number of water columns
nForcPeriods = length(forc.month); % number of distinct periods of forcing data
months = forc.month;
nvars = init.nvars;
init = rmfield(init, 'nvars');
varNames = fieldnames(init);
nEquations = nvars * ndepth; % number of (discretised) equations to solve
isLand = domain.isLand;

isLeapYear = false;
dv = datevec(forc.time);
isFeb = dv(:,2) == 2;
if any(isFeb)
    yr = dv(isFeb,1);
    dv = datevec(datenum([yr, 2, 29]));
    if dv(2) == 2, isLeapYear = true; end
end

% For efficient use in the parfor loop, reconfigure structure of initial
% values by permuting unique water columns to the trailing dimension
init = structfun(@(z) permute(z, [3, 1, 2]), init, 'UniformOutput', false);
init = structfun(@(z) z(:,:), init, 'UniformOutput', false);
init = cell2mat(struct2cell(init));
% yy = reshape(permute(init_, [2, 1]), [nColumns, domain.ndepth, nvars]);
% yy = reshape(permute(yy, [3, 2, 1]), nvars, domain.ndepth, domain.nlon, domain.nlat);
% yy = permute(yy, [1, 3, 4, 2]);


% Reconfigure structure of forcing data for use in the parfor loop
forc_ = forc;
fields = fieldnames(forc);
forcDims = structfun(@(z) ndims(z), forc);
forc = rmfield(forc, fields(~ismember(forcDims, [3, 4])));
forc = structfun(@(z) permute(z, [3, 4, 1, 2]), forc, 'UniformOutput', false);
forc = structfun(@(z) z(:,:,:), forc, 'UniformOutput', false);

% Pre-allocate output array -- construction styled for compatibiliy with
% parallel loop
daysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31];
if isLeapYear
    daysPerMonth(2) = 29;
end
secsPerMonth = 86400 * daysPerMonth; % length of each month
secsPerMonth = secsPerMonth(months); % retain modelled months only
nOutPerMonth = secsPerMonth / pars.dt_out; % number of model outputs each month
nOutputTimes = sum(nOutPerMonth);
xx = cell(nForcPeriods, 1);
for k = 1:nForcPeriods
    xx{k} = nan(nEquations, nOutPerMonth(k));
end
OUT = repmat(xx, 1, nColumns);


%% Integration

% Loop through valid horizontal grid cells

% Index valid horizontal grid cells
ind = domain.horizGridIndex(:,:,1) .* ~isLand;
ind(ind == 0) = [];
nvalid = length(ind);
OUT_ = OUT(:,ind);
init_ = init(:,ind);
% Loop through each water column
parfor i = 1:nvalid
    ii = ind(i);
    % Forcing data
    Forcing = structfun(@(z) z(:,:,ii), forc, 'UniformOutput', false);
    % Initial state
    y0 = init_(:,i);
    % Integrating method
    odeSolve = str2func(odeIntegrator);
    % Integration times -- start/end times matching months
    T0 = zeros(1, nForcPeriods); % integration start times
    T1 = secsPerMonth; % and end times (seconds)
    % Integrate step-wise between successive periods of forcing data
    for j = 1:nForcPeriods
        t0 = T0(j);
        t1 = T1(j);
        forcing = structfun(@(z) z(:,j), Forcing, 'UniformOutput', false);
        % Integrate
        sol = odeSolve(@(t, y) risk_map_model(t, y, domain, pars, forcing), ...
            [t0 t1], y0, odeOptions);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % This method treats forcing data as constant within each month. It
        % could be interpolated across months by including time-step j as
        % an argument to risk_map_model.
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Store solutions for selected output times
        outTimes = pars.dt_out:pars.dt_out:t1;
        OUT_{j,i} = deval(sol, outTimes);
        % Update initials for next integration increment
        y0 = OUT_{j,i}(:,end);
    end
end

%% Tidy up

% Move solutions into full output array
for i = 1:nvalid
    ii = ind(i);
    for j = 1:nForcPeriods
        OUT{j,ii} = OUT_{j,i};
    end
end

% Convert solution stored in cell-array OUT into more readable format
x = OUT(:);
x = cell2mat(x');
x = reshape(x, nEquations, nOutputTimes, nColumns);
x = permute(x, [2, 3, 1]);
x = reshape(x, nOutputTimes, nColumns, ndepth, nvars);
x = permute(x, [4, 3, 1, 2]);
x = reshape(x, [nvars, ndepth, nOutputTimes, nlon, nlat]);
x = permute(x, [1, 4, 5, 2, 3]);
% Output array [variable, lon, lat, depth, time] contains complete solution
out = x;

% Optional output struct stores solution in more readable format
outStruct.lon = domain.lon;
outStruct.lat = domain.lat;
outStruct.depth = domain.depth;
% Include output times
initTime = datevec(forc_.time(1));
initTime(3) = 1;
initTime(4:6) = 0;
timeIncrements = linspace(pars.dt_out, pars.dt_out * nOutputTimes, nOutputTimes);
timeIncrements = [zeros(nOutputTimes, 5) timeIncrements'];
allTimes = initTime + timeIncrements;
allTimes = datenum(allTimes);
% datestr(datenum(allTimes))
outStruct.time = allTimes;
x = permute(x, [2, 3, 4, 5, 1]);
for i = 1:nvars
    outStruct.(varNames{i}) = x(:,:,:,:,i);
end

