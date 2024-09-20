function load_constants_2D()

% I have used the same constants as Mirza
n0      = 1;                     % Refractive index of air
ns      = 0.571 - 1i*0.636;      % Refractive index of the dielectric substrate
n1      = 5.05506 - 1i*3.20418;  % Refractive index of the protein
eta0    = 377;                   % 377 Ohm, impedance of free space
lambda0 = 325*10^(-9);           % Wavelength in free space
Gamma_r = (n0 - ns)/(n0 + ns);   % Fresnel reflection coefficient
Gamma_t = 2*n0/(n0 + ns);        % Fresnel transmission coefficient
k0      = 2*pi/lambda0;          % Wave number in free space
ks      = k0*ns/n0;              % Wave number in substrate
N       = 100;                   % Number of auxiliary sources, either interior or exterior, or test points per segment
%nw.alpha = 0.86;                 % Scaling for placement of auxiliary sources