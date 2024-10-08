function Hx_tot_inside = Hx_tot_inside_matrix(coord,coord_ext)

% Get coordinates
x     = coord.x;
y     = coord.y;
x_ext = coord_ext.x;
y_ext = coord_ext.y;

% Load constants
const   = load_constants();
k1      = const.k1;
n1      = const.n1;
eta0    = const.eta0;

% Do precomputations
H12 = @(z) besselh(1,2,z);
numerical = @(x,y) sqrt(x.^2+y.^2);
abs_ext     = numerical(x - x_ext, y - y_ext);

% Compute field
Hx_tot_inside = - 1i*n1/eta0 * 1./abs_ext .* (y - y_ext) .* H12(k1*abs_ext);