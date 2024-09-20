function Ez_tot_inside = Ez_tot_inside_matrix(x,y,x_ext,y_ext)
% x and y has to be vectors of same length
assert(length(x)==length(y),"x and y does not have the same length!")
assert(length(x_ext)==length(y_ext),"x' and y' does not have the same length!")

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants_nanowires();

% Do precomputations
H02 = @(z) besselh(0,2,z);
numerical = @(x,y) sqrt(x.^2+y.^2);
abs_ext     = numerical(x - x_ext, y - y_ext);

% Compute field
Ez_tot_inside = H02(k1*abs_ext);