function [EHx, EHy, EHz] = EH(x,y,z)

% Get constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, k1] = load_constants();

% What should these numbers be?
I = 1;
dl = 1;
k = 1;

% Compute variables in the spherical coordinate system with origo at (x_int,y_int,z_int)
rho = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z./rho);
phi = pi + sign(y) .* (acos(x./rho) - pi);

% Compute fields in the spherical coordinate system
EH_rho   = 1i * eta0 * I * dl/(4*pi) * (k./rho - 1i./rho.^2 - 1./(k*rho.^3)).*exp(-1i*k*rho).*sin(theta);
EH_theta =      eta0 * I * dl/(4*pi) * (1./rho.^2 - 1i./(k*rho.^3)).*exp(-1i*k*rho).*cos(theta);

% Compute fields
EHx = EH_rho .* sin(theta) .* cos(phi) + EH_theta .* cos(theta) .* cos(phi);
EHy = EH_rho .* sin(theta) .* sin(phi) + EH_theta .* cos(theta) .* sin(phi);
EHz = EH_rho .* cos(theta)             - EH_theta .* sin(theta);