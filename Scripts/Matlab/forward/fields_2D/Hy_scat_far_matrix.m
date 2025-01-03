function Hy_scat = Hy_scat_far_matrix(coord,coord_int)

% Get coordinates
x     = coord.x;
y     = coord.y;
x_int = coord_int.x;
y_int = coord_int.y;

% Load constants
const   = load_constants();
Gamma_r = const.Gamma_r;
k0      = const.k0;
eta0    = const.eta0;

% Do precomputations
alpha = 1; 
H12 = @(z) sqrt(2./(pi*z)) .* exp(-1i * (z - alpha*pi/2 - pi/4)); 
numerical = @(x,y) sqrt(x.^2+y.^2);
abs_int     = numerical(x - x_int, y - y_int);
abs_int_ref = numerical(x - x_int, y + y_int);

% Compute field
Hy_scat = 1i/eta0 * (x - x_int) .* (1./abs_int     .* H12(k0*abs_int) + ...
                          Gamma_r * 1./abs_int_ref .* H12(k0*abs_int_ref));