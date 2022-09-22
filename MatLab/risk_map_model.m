function dvdt = risk_map_model(t, y, domain, pars, forc)
% Differential equations describing flux of faecal pellet carbon.

%% Production
p = production(Z,P);

%% Remineralisation
r = remineralisation(M);

%% Sinking
s = sinking(M);

%% Total flux
dvdt = p - r - s;


end


function out = production(z, p)

end

function out = remineralisation(m)

end

function out = sinking(m)

end