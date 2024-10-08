function H_ref = H_ref_vector(x,y,z)

H_ref = zeros(length(x),3);

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, k1] = load_constants();

% Compute field
H_ref(:,1) = 1/eta0 * Gamma_r * exp(-1i*k0*y);