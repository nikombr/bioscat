function [Etot_inside, Htot_inside] = interior_fields(coord, coord_ext, scenario, lambda0, far_field_approximation)

% Get constants
const    = load_constants(lambda0);
k0       = const.k0;
epsilon0 = const.epsilon0;
mu0      = const.mu0;

% Get coordinates
x     = coord.x;
y     = coord.y;
x_ext = coord_ext.x;
y_ext = coord_ext.y;

% Allocation
m           = length(x);
n           = length(x_int);
Etot_inside = zeros(m,n,3);
Htot_inside = zeros(m,n,3);

% Do precomputations
if far_field_approximation
    alpha = 0; 
    H02 = @(z) sqrt(2./(pi*z)) .* exp(-1i * (z - alpha*pi/2 - pi/4)); 
    alpha = 1; 
    H12 = @(z) sqrt(2./(pi*z)) .* exp(-1i * (z - alpha*pi/2 - pi/4));
else
    H02 = @(z) besselh(0,2,z); % Hankel function of zero order and of the second kind
    H12 = @(z) besselh(1,2,z); % Hankel function of first order and of the second kind
end

numerical = @(x,y) sqrt(x.^2 + y.^2);
abs_ext     = numerical(x - x_ext, y - y_ext);

if scenario == 1

elseif scenario == 2 


end