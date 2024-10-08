function Ez_scat = Ez_scat_matrix(coord,coord_int)

% Get coordinates
x     = coord.x;
y     = coord.y;
x_int = coord_int.x;
y_int = coord_int.y;

% Load constants
const   = load_constants();
Gamma_r = const.Gamma_r;
k0      = const.k0;

% Do precomputations
H02 = @(z) besselh(0,2,z);
numerical = @(x,y) sqrt(x.^2+y.^2);
abs_int     = numerical(x - x_int, y - y_int);
abs_int_ref = numerical(x - x_int, y + y_int);

% Compute field
Ez_scat = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);