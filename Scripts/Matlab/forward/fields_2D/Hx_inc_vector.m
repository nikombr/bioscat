function Hx_inc = Hx_inc_vector(x,y)
% x and y has to be vectors of same length
assert(length(x)==length(y),"x and y does not have the same length!")

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants();

% Compute field
Hx_inc = - 1/eta0 * exp(1i*k0*y);