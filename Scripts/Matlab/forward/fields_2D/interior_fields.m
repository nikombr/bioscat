function [Etot_inside, Htot_inside] = interior_fields(coord, coord_ext, scenario, lambda0)

% Get constants
const    = load_constants(lambda0);
k1       = const.k1;
n1       = const.n1;
epsilon0 = const.epsilon0;
mu0      = const.mu0;
eta0     = const.eta0;

% Get coordinates
x     = coord.x;
y     = coord.y;
x_ext = coord_ext.x;
y_ext = coord_ext.y;

% Allocation
m           = length(x);
n           = length(x_ext);
Etot_inside = zeros(m,n,3);
Htot_inside = zeros(m,n,3);

% Do precomputations
H02 = @(z) besselh(0,2,z); % Hankel function of zero order and of the second kind
H12 = @(z) besselh(1,2,z); % Hankel function of first order and of the second kind

numerical = @(x,y) sqrt(x.^2 + y.^2);
abs_ext     = numerical(x - x_ext, y - y_ext);

% Do precomputations of functions used twice
H12_eval     = H12(k1*abs_ext);

if scenario == 1

    Etot_inside(:,:,3) = H02(k1*abs_ext);

    Htot_inside(:,:,1) = - 1i*n1/eta0 * 1./abs_ext .* (y - y_ext) .* H12_eval;

    Htot_inside(:,:,2) =   1i*n1/eta0 * 1./abs_ext .* (x - x_ext) .* H12_eval;

elseif scenario == 2 

    Etot_inside(:,:,1) =   1i*n1*mu0/(eta0*epsilon0) * 1./abs_ext .* (y - y_ext) .* H12_eval;

    Etot_inside(:,:,2) = - 1i*n1*mu0/(eta0*epsilon0) * 1./abs_ext .* (x - x_ext) .* H12_eval;

    Htot_inside(:,:,3) = H02(k1*abs_ext);

end