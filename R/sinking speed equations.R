# Sinking speed equations

ss_cylinder = function(L, D, rho, rho_p, mu, g = 981){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Terminal sinking speed of cylinders in Stokes' flow regime.
  # Equation uses gram-centimetre-second units.
  # See Komar et al., 1981 -- doi:10.4319/lo.1981.26.1.0172
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # particle length: [L] = cm
  # particle diameter: [D] = cm
  # seawater density: [rho] = g / cm^3 (= 1e3 kg / m^3)
  # particle density: [rho_p] = g / cm^3
  # seawater viscosity: [mu] = g / cm / s
  # acceleration due to gravity: [g] = cm / s^2
  # returns sink speed: [ss] = cm / s
  0.079 / mu * {rho_p - rho} * g * L^2 * {L / D}^-1.664
}

shape_ellipsoid = function(Ds, Di, Dl){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Dimensionless shape parameter for ellipsoidal particles.
  # Used in sinking equation.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # particle short, medium, long diameters: [Ds,Di,Dl] = cm
  Ds * {{{Ds ^ 2 + Di ^ 2 + Dl ^2} / 3} ^ -0.5}
}

ss_ellipsoid = function(Dn, Ds, Di, Dl, rho, rho_p, mu, g = 981){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Terminal sinking speed of ellipsoids in Stokes' flow regime.
  # Equation uses gram-centimetre-second units.
  # See Komar et al., 1981 -- doi:10.4319/lo.1981.26.1.0172
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # particle nominal diameter: [Dn] = cm
  # particle short, medium, long diameters: [Ds,Di,Dl] = cm
  # seawater density: [rho] = g / cm^3 (= 1e3 kg / m^3)
  # particle density: [rho_p] = g / cm^3
  # seawater viscosity: [mu] = g / cm / s
  # acceleration due to gravity: [g] = cm / s^2
  # returns sink speed: [ss] = cm / s
  E = shape_ellipsoid(Ds, Di, Dl)
  1 / 18 / mu * {rho_p - rho} * g * Dn ^ 2 * E ^ 0.38
}





