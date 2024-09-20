function nw = setup_nanowire(xc, r)

nw = struct;

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants_nanowires();

if nargin == 0
    % Default values if we have no input
    nw.r = 0.5*lambda0; % Radius of nanowire
    nw.xc = 0;          % Placement of nanowire
end

% Set default values for variables used later on
nw.C = [];                      % Complex amplititudes for computing the exterior field
nw.D = [];                      % Complex amplititudes for computing the interior field
nw.x = []; nw.y = [];           % Test points on the boundary
nw.x_int = []; nw.y_int = [];   % Auxiliary sources in the interior
nw.x_ext = []; nw.y_ext = [];   % Auxiliary sources in the exterior