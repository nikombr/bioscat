function nw = setup_nanowire(N, xc, r)

nw = struct;

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants();

if nargin < 2
    % Default values if we have no input
    r  = 0.5*lambda0; % Radius of nanowire;
    xc = 0;          % Placement of nanowire
end
nw.r  = r; % Radius of nanowire
nw.xc = xc;          % Placement of nanowire

% Compute auxiliary sources and test points
nw.phi = linspace(0,2*pi,N)';
nw.x = r*cos(nw.phi) + xc;
nw.y = r*sin(nw.phi) + r;
nw.x_int = alpha*r*cos(nw.phi') + xc;
nw.y_int = alpha*r*sin(nw.phi') + r;
nw.x_ext = 1/alpha*r*cos(nw.phi') + xc;
nw.y_ext = 1/alpha*r*sin(nw.phi') + r;

% Set default values for variables used later on
nw.C = [];                      % Complex amplititudes for computing the exterior field
nw.D = [];                      % Complex amplititudes for computing the interior field

