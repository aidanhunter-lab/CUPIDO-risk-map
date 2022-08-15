function v = sinking_speed(shape, rho, rho_p, mu, g, varargin)
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
% varargin contains faecal pellet dimensions -- these depend on shape [cm]
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extractVarargin(varargin)

switch shape

    case 'cylinder'
        % Pellet length and diameter must be given
        if ~exist('L', 'var') || ~exist('D', 'var')
            error('Faecal pellet length and diameter, L and D, must be specified in varargin.')
        end
        v = 0.079 ./ mu .* (rho_p - rho) .* g .* L ^ 2 * (L ./ D) ^ -1.664;

    case 'ellipsoid'
        % Pellet small, medium and large diameters must be given
        if ~exist('Ds', 'var') || ~exist('Di', 'var') || ~exist('Dl', 'var')
            error('Faecal pellet diameters (short = Ds, intermediate = Di, large = Dl) must be specified in varargin.')
        end
        E = shape_ellipsoid(Ds, Di, Dl);
        v = 1 ./ 18 ./ mu .* {rho_p - rho} .* g .* Dn ^ 2 .* E ^ 0.38;
end
end

function s = shape_ellipsoid(Ds, Di, Dl)
% Dimensionless shape parameter for ellipsoidal particles.
% Particle short, medium, long diameters: [Ds,Di,Dl] = cm
s = Ds .* (((Ds ^ 2 + Di ^ 2 + Dl ^ 2) ./ 3) ^ -0.5);
end
