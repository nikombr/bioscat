function const = load_constants()

const = struct;

% I have used the same constants as Mirza
const.ns      = 1.6; %0.571 - 1i*0.636;% Refractive index of the dielectric substrate
const.eta0    = 377;                   % 377 Ohm, impedance of free space
const.n0      = 1;                     % Refractive index of air
const.lambda0 = 325*10^(-9);           % Wavelength in free space
const.Gamma_r = (const.n0 - const.ns)/(const.n0 + const.ns);   % Fresnel reflection coefficient
const.Gamma_t = 2*const.n0/(const.n0 + const.ns);        % Fresnel transmission coefficient
const.k0      = 2*pi/const.lambda0;          % Wave number in free space
const.ks      = const.k0*const.ns/const.n0;              % Wave number in substrate

% We have the same material and computationally properties for all the nanowires
const.n1    = 1.6; %5.05506 - 1i*3.20418; % Refractive index of the nanowire
const.alpha = 0.86;  % Scaling for placement of auxiliary sources

% Compute constants from the other
const.k1 = const.k0*const.n1/const.n0; % Wave number in nanowire