function const = load_constants(lambda0)

const = struct;

% Physical constants
const.eta0     = 377;                      % 377 Ohm, impedance of free space
const.n0       = 1;                        % Refractive index of air
const.epsilon0 = 8.8541878188  * 10^(-12); % Permittivity of free space
const.mu0      = 1.25663706127 * 10^(-6);  % Permeability of free space

% Chosen constants (I have used the same constants as Mirza for some of them)
if nargin == 0
    const.lambda0 = 325*10^(-9); % Default value of wavelength in free space
else
    const.lambda0 = lambda0;     % Chosen value of wavelength in free space
end
const.ns    = 1.6;  % Refractive index of the substrate %0.571 - 1i*0.636;
const.n1    = 1.6;  % Refractive index of the nanowire/nanostructure % 5.05506 - 1i*3.20418;
const.alpha = 0.86; % Scaling for placement of auxiliary sources for the nanowires

% Computed constants
const.Gamma_r = (const.n0 - const.ns)/(const.n0 + const.ns); % Fresnel reflection coefficient
const.Gamma_t = 2*const.n0/(const.n0 + const.ns);            % Fresnel transmission coefficient
const.k0      = 2*pi/const.lambda0;                          % Wave number in free space
const.ks      = const.k0*const.ns/const.n0;                  % Wave number in substrate
const.k1      = const.k0*const.n1/const.n0;                  % Wave number in nanowire/nanostructure

