function Ez_scat = Ez_scat_matrix(x,y,x_int,y_int)
% x and y has to be vectors of same length
assert(length(x)==length(y),"x and y does not have the same length!")
assert(length(x_int)==length(y_int),"x' and y' does not have the same length!")

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants();

% Do precomputations
H02 = @(z) besselh(0,2,z);
numerical = @(x,y) sqrt(x.^2+y.^2);
abs_int     = numerical(x - x_int, y - y_int);
abs_int_ref = numerical(x - x_int, y + y_int);

% Compute field
Ez_scat = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);