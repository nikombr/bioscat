function [Escat, Hscat] = scattered_fields(coord, coord_int, scenario, lambda0, far_field_approximation)

% Get constants
const    = load_constants(lambda0);
k0       = const.k0;
epsilon0 = const.epsilon0;
mu0      = const.mu0;

% Get coordinates
x     = coord.x;
y     = coord.y;
x_int = coord_int.x;
y_int = coord_int.y;

% Allocation
m = length(x);
n = length(x_int);
Escat = zeros(m,n,3);
Hscat = zeros(m,n,3);

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
abs_int     = numerical(x - x_int, y - y_int);
abs_int_ref = numerical(x - x_int, y + y_int);

% Do precomputations of functions used twice
H12_eval     = H12(k0*abs_int);
H12_eval_ref = H12(k0*abs_int_ref);

if scenario == 1
    
    % z-component of the E-field
    Escat(:,:,3) = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);
    
    % x-component of the H-field
    Hscat(:,:,1) = - 1i/eta0 * (1./abs_int     .* H12_eval      .* (y - y_int) + ...
                      Gamma_r * 1./abs_int_ref .* H12_eval_ref  .* (y + y_int));
    
    % y-component of the H-field
    Hscat(:,:,2) = 1i/eta0 * (x - x_int) .* (1./abs_int     .* H12_eval + ...
                                   Gamma_r * 1./abs_int_ref .* H12_eval_ref);

elseif scenario == 2 

    % x-component of the H-field
    Escat(:,:,1) = 1i*mu0/(eta0*epsilon0) * (1./abs_int     .* H12_eval      .* (y - y_int) + ...
                                   Gamma_r * 1./abs_int_ref .* H12_eval_ref  .* (y + y_int));
    
    % y-component of the H-field
    Escat(:,:,2) = - 1i*mu0/(eta0*epsilon0) * (x - x_int) .* (1./abs_int     .* H12_eval + ...
                                                    Gamma_r * 1./abs_int_ref .* H12_eval_ref);

    % z-component of the E-field
    Hscat(:,:,3) = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);

end