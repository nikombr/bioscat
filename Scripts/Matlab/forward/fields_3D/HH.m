function [HHx, HHy, HHz] = HH(x,y,z)

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
HH_phi = 1i * I * dl/(4*pi) * (k./rho - 1i./rho.^2) .* exp(-1i*k*rho) .* sin(theta);

% Compute fields
HHx = - HH_phi .* sin(phi);
HHy =   HH_phi .* cos(phi);
HHz =   HH_phi * 0;