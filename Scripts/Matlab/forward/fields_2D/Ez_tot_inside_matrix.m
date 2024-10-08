function Ez_tot_inside = Ez_tot_inside_matrix(coord,coord_ext)

% Get coordinates
x     = coord.x;
y     = coord.y;
x_ext = coord_ext.x;
y_ext = coord_ext.y;

% Load constants
const   = load_constants();
k1      = const.k1;

% Do precomputations
H02         = @(z) besselh(0,2,z);
numerical   = @(x,y) sqrt(x.^2+y.^2);
abs_ext     = numerical(x - x_ext, y - y_ext);

% Compute field
Ez_tot_inside = H02(k1*abs_ext);