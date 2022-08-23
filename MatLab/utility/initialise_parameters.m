function pars = initialise_parameters()

pars.mu = 0.00189 * 1e1;  % seawater viscosity, g / cm / s. [1 N s / m^2 = 1 Pa s = 10 g / cm / s]
pars.g = 9.81 * 1e2;      % acceleration due to gravity, cm / s^2. [1 m / s^2 = 10^2 cm / s^2]

% See Atkinson et al., 2012, for krill faecal pellet properties
pars.rho_p = 1116 * 1e-3; % faecal pellet density, g / cm^3. [1 kg / m^3 = 10^-3 g / cm^3]
% pars.rho_p = 1220 * 1e-3; % value from Komar 1980 

% For now, let's assume that faecal pellets are either cylindrical OR
% ellipsoidal. We may introduce a mixture of shapes later...
pars.shape = 'cylinder';  % faecal pellet shape (may be cylinder or ellipsoid)

% It seems that Atkinson's data refer to cylindrical pellets
pars.L = 2667 * 1e-4; % median length (cm)
pars.D = 178 * 1e-4; % median width (diameter) (cm)

switch pars.shape
    case 'cylinder'
        pars.Ds = []; pars.Di = []; pars.Dl = [];
        pars.V = 0.25 * pi * pars.D ^ 2 * pars.L;
    case 'ellipsoid'
        % assume that the two smallest diameters are equal
        pars.Ds = pars.D; pars.Di = pars.D;
        pars.Dl = pars.L;
        pars.V = 1 / 6 * pi * pars.Ds * pars.Di * pars.Dl;
end

% The faecal pellet volumes from the Atkinson data are at the uppper end of
% the range given for euphausiid pellets (approx 1e-5 -> 5e-5) in 
% Komar et al (1981). This makes sense as Antarctic krill are largest
% euphausiid. The volume is also close to the standard volume (5e-5) used
% by Schmidt et al (2012).
