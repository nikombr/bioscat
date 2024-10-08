function H_inc = H_inc_vector(x,y,z)

H_inc = zeros(length(x),3);

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, k1] = load_constants();

% Compute field
H_inc(:,1) = - 1/eta0 * exp(1i*k0*y);