function [v, Re] = sinking_speed(shape, rho, rho_p, mu, g, varargin)
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
if ~exist('useMeans_Re', 'var'), useMeans_Re = false; end
if ~exist('returnMax_Re', 'var'), returnMax_Re = false; end

switch shape

    case 'cylinder'
        % Pellet length and diameter must be given
        if isempty(L) || isempty(D)
            error('Faecal pellet length and diameter, L and D, must be specified in varargin.')
        end

        v = sink_speed_cylinder(rho_p, rho, mu, g, L, D);

    case 'ellipsoid'
        % Pellet small, medium and large diameters must be given
        if isempty(Ds) || isempty(Di) || isempty(Dl)
            error('Faecal pellet diameters (short = Ds, intermediate = Di, large = Dl) must be specified in varargin.')
        end

        [v, Dn] = sink_speed_ellipsoid(rho_p, rho, mu, g, Ds, Di, Dl);

end

% Calculate Re only if it is requested as output argument

if nargout > 1
    switch shape, case 'cylinder'
        Dn = nominal_diameter_cylinder(L, D);
    end
    if useMeans_Re && returnMax_Re
        Re = Reynolds_number(rho, mu, v, Dn, 'max', 'mean');
    elseif useMeans_Re && ~returnMax_Re
        Re = Reynolds_number(rho, mu, v, Dn, 'mean');
    elseif ~useMeans_Re && returnMax_Re
        Re = Reynolds_number(rho, mu, v, Dn, 'max');
    elseif ~useMeans_Re && ~returnMax_Re
        Re = Reynolds_number(rho, mu, v, Dn);
    end
end

end % END OF MAIN FUNCTION

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function s = shape_ellipsoid(Ds, Di, Dl)
% Dimensionless shape parameter for ellipsoidal particles.
% Particle short, medium, long diameters: [Ds,Di,Dl] = cm
s = Ds .* (((Ds ^ 2 + Di ^ 2 + Dl ^ 2) ./ 3) ^ -0.5);
end

function Dn = nominal_diameter_cylinder(L, D)
% Diameter of sphere with volume equal to volume of cylinder described by
% length, L, and diameter, D.
V = 0.25 * pi * D ^ 2 * L; % cylinder volume
Dn = (6 * V / pi) ^ (1/3);
end

function Dn = nominal_diameter_ellipsoid(Ds, Di, Dl)
% Diameter of sphere with volume equal to volume of ellipsoid described by
% principal axes Ds, Di and Dl.
V = pi / 6 .* Ds .* Di .* Dl; % ellipsoid volume
Dn = (6 * V / pi) ^ (1/3);
end

function v = sink_speed_cylinder(rho_p, rho, mu, g, L, D)
v = 0.079 ./ mu .* (rho_p - rho) .* g .* L ^ 2 * (L ./ D) ^ -1.664;
end

function [v, Dn] = sink_speed_ellipsoid(rho_p, rho, mu, g, Ds, Di, Dl)
E = shape_ellipsoid(Ds, Di, Dl);
Dn = nominal_diameter_ellipsoid(Ds, Di, Dl);
v = 1 ./ 18 ./ mu .* (rho_p - rho) .* g .* Dn ^ 2 .* E ^ 0.38;
end

function Re = Reynolds_number(rho, mu, v, Dn, varargin)
% Estimate the maximum Reynolds number to test for laminar flow regime.
% Optional arguments may include 'mean' and/or 'max'. It is only the
% maximum Reynolds number that is required to test the assumption of
% laminar flow. Taking averages reduce RAM requirements -- Re will vary
% most with depth so average over lat-lon-month.
useMeans = any(strcmp('mean', varargin));
returnMax = any(strcmp('max', varargin));
if useMeans
    rho = mean(rho, [1, 2, 4], 'omitnan');
    v = mean(v, [1, 2, 4], 'omitnan');
end
Re = rho .* v .* Dn ./ mu;
if returnMax
    Re = max(Re(:));
end
end

