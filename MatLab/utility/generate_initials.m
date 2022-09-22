function v = generate_initials(domain, pars, forc)
% Create initial values for all variables in differential equations.

% Faecal pellet production is some function of krill abundance and food
% availability, so use these data to derive initials...
% forc.krill

% For now, just to move forward quickly with model building, assume that
% initial state is a single faecal pellet per grid cell.

% WHY IS THE BOTTOM LAYER EMPTY??? THIS SUGGESTS SOMETHING'S AMISS WITH THE
% DOMAIN FUNCTION...
nFP = ones(domain.gridSize) .* domain.isWaterColumn; % number of faecal pellets

VFP = pars.V .* nFP; % total volume of faecal pellets (cm^3)
CFP = pars.fpCm_summer .* VFP; % total carbon mass of faecal pellets (mg C)

% Store initial values for all variables
v_.CFP = CFP;

fields = fieldnames(v_);
v.nvars = length(fields);
for i = 1:v.nvars
    v.(fields{i}) = v_.(fields{i});
end
