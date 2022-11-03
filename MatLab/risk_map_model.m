function dvdt = risk_map_model(~, y, domain, pars, forcing)
% Differential equations describing flux of faecal pellet carbon.

%% Initial conditions
CFP = y; % Total faecal pellet carbon [g C]

%% Production
% p = production(Z,P);

% A = 0.5 .* (domain.area(1:end-1) + domain.area(2:end)); % average cross-sectional area [m^2] of grid cell

krill = forcing.krill .* domain.volume; % total krill [number of individuals] in grid cell -- should be a vector over depth, where only the uppermost layer is populated


% Convert number of individuals into faecal pellet production. I could use
% some assumption of average size of krill, or use the length data. It's
% probably best to use the length data, but for starters just make some
% assumption about averge size...
% This will involve converting krill abundance from from numbers per m2
% into numbers per grid cell using the domain cell areas

% faecal pellet production if a function of krill biomass...

krillDryWtMean = 10^-3 .* krill .* exp(forcing.krillDryWtMeanLog); % total krill [g dry mass] in grid cell

FPprod = pars.FPprod; % mm^3 / h / g dry mass
FPprod = 10^-3 ./ 3600 .* FPprod; % cm^3 / s / g dry mass
FPprod = pars.fpCm_summer .* FPprod; % g C / s / g dry mass

p = FPprod .* krillDryWtMean; % g C / s
p(isnan(p)) = 0;


%% Remineralisation
% r = remineralisation(M);

%% Sinking
% s = sinking(M);
w = sink_speed(pars.shape, 1e-3 * forcing.density_seawater, pars.rho_f, ...
    pars.mu, pars.g, 'L', pars.L, 'D', pars.D); % [cm / s]
w = 10^-2 .* w; % [m / s]

% s = sink_flux(CFP, w, diff(100 * domain.depthgrid));
s = sink_flux(CFP, w, diff(domain.depthgrid)); % [g C / s]


%% Total flux
% dvdt = p - r - s;
dvdt = p + s;

end

%% Production equations
% function out = production(z, p)
% end

%% Remineralisation equations
% function out = remineralisation(m)
% end



%% Sinking equations
% function v = sink_speed(m)
% end

function [v, Re] = sink_speed(shape, rho, rho_p, mu, g, varargin)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Terminal sinking speed equations for krill faecal pellets in Stokes flow
% regime. Equations use gram-centimetre-second units.
% See Komar et al., 1981 -- doi:10.4319/lo.1981.26.1.0172
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input arguments:
% shape my be 'cylinder' or ellipsoid.
% rho = seawater density [g / cm^3 (= 1e3 kg / m^3)]
% rho_p = particle density [g / cm^3]
% mu = seawater viscosity [g / cm / s]
% g = acceleration due to gravity [cm / s^2]
% varargin must contain faecal pellet dimensions [cm] (dependent on shape),
% and may optionally contain arguments augmenting the Re calculation for
% numerical efficiency and reducing RAM requirements.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extractVarargin(varargin)
if ~exist('L', 'var'), L = []; end
if ~exist('D', 'var'), D = []; end
if ~exist('Ds', 'var'), Ds = []; end
if ~exist('Di', 'var'), Di = []; end
if ~exist('Dl', 'var'), Dl = []; end
% if ~exist('useMeans_Re', 'var'), useMeans_Re = false; end
if ~exist('returnMax_Re', 'var'), returnMax_Re = false; end
switch shape
    case 'cylinder'
        if isempty(L) || isempty(D) % Pellet length and diameter must be given
            error('Faecal pellet length and diameter, L and D, must be specified in varargin.')
        end
        v = sink_speed_cylinder;
    case 'ellipsoid'
        if isempty(Ds) || isempty(Di) || isempty(Dl) % Pellet small, medium and large diameters must be given
            error('Faecal pellet diameters (short = Ds, intermediate = Di, large = Dl) must be specified in varargin.')
        end
        [v, Dn] = sink_speed_ellipsoid;
end
% Calculate Re only if it is requested as output argument
if nargout > 1
    switch shape, case 'cylinder'
        Dn = nominal_diameter_cylinder;
    end
    if ~returnMax_Re
        Re = Reynolds_number;
    else
        Re = Reynolds_number('max');
    end
end
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% END OF MAIN FUNCTION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nested functions
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function V = volume_cylinder, V = 0.25 * pi * D ^ 2 * L; end
    function V = volume_ellipsoid, V = pi / 6 .* Ds .* Di .* Dl; end
    function s = shape_ellipsoid
        % Dimensionless shape parameter for ellipsoidal particles.
        % Particle short, medium, long diameters: [Ds,Di,Dl] = cm
        s = Ds .* (((Ds ^ 2 + Di ^ 2 + Dl ^ 2) ./ 3) ^ -0.5);
    end
    function Dn = nominal_diameter_cylinder
        % Diameter of sphere with volume equal to volume of cylinder described by
        % length, L, and diameter, D.
        V = volume_cylinder;
        Dn = (6 * V / pi) ^ (1/3);
    end
    function Dn = nominal_diameter_ellipsoid
        % Diameter of sphere with volume equal to volume of ellipsoid described by
        % principal axes Ds, Di and Dl.
        V = volume_ellipsoid;
        Dn = (6 * V / pi) ^ (1/3);
    end
    function v = sink_speed_cylinder
        v = 0.079 ./ mu .* (rho_p - rho) .* g .* L ^ 2 * (L ./ D) ^ -1.664;
    end
    function [v, Dn] = sink_speed_ellipsoid
        E = shape_ellipsoid;
        Dn = nominal_diameter_ellipsoid;
        v = 1 ./ 18 ./ mu .* (rho_p - rho) .* g .* Dn ^ 2 .* E ^ 0.38;
    end
    function Re = Reynolds_number(varargin)
        % Estimate the maximum Reynolds number to test for laminar flow regime.
        returnMax = any(strcmp('max', varargin));
        Re = rho .* v .* Dn ./ mu;
        if returnMax
            Re = max(Re(:));
        end
    end
end

function v = sink_flux(u,s,w)
% Rate of change of u due to sinking.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         s = sinking speed, size(s)=[1 nvar]
%         w = depth layer widths, size(w)=[nz 1]
v = ((-s) ./ w) .* diff([zeros(1,size(u,2)); u]);
end



