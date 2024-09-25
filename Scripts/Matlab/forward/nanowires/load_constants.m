function [eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants()

% I have used the same constants as Mirza
ns      = 1.6; %0.571 - 1i*0.636;% Refractive index of the dielectric substrate
eta0    = 377;                   % 377 Ohm, impedance of free space
n0      = 1;                     % Refractive index of air
lambda0 = 325*10^(-9);           % Wavelength in free space
Gamma_r = (n0 - ns)/(n0 + ns);   % Fresnel reflection coefficient
Gamma_t = 2*n0/(n0 + ns);        % Fresnel transmission coefficient
k0      = 2*pi/lambda0;          % Wave number in free space
ks      = k0*ns/n0;              % Wave number in substrate

% We have the same material and computationally properties for all the nanowires
n1    = 1.6; %5.05506 - 1i*3.20418; % Refractive index of the nanowire
alpha = 0.86;                 % Scaling for placement of auxiliary sources

% Compute constants from the other
k1 = k0*n1/n0; % Wave number in nanowire