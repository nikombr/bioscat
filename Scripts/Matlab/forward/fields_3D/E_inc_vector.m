function E_inc = E_inc_vector(x,y,z)

E_inc = zeros(length(x),3);

[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, k1] = load_constants();

% Compute field
E_inc(:,3) = exp(1i*k0*y);