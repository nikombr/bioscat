function Hy_scat = Hy_scat_matrix(x,y,x_int,y_int)
% x and y has to be vectors of same length
assert(length(x)==length(y),"x and y does not have the same length!")
assert(length(x_int)==length(y_int),"x' and y' does not have the same length!")

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants_nanowires();

% Do precomputations
H12 = @(z) besselh(1,2,z); 
numerical = @(x,y) sqrt(x.^2+y.^2);
abs_int     = numerical(x - x_int, y - y_int);
abs_int_ref = numerical(x - x_int, y + y_int);

% Compute field
Hy_scat = - 1i/eta0 * (x - x_int) .* (1./abs_int     .* H12(k0*abs_int) + ...
                                      1./abs_int_ref .* H12(k0*abs_int_ref));