function L = nanostructures2Dsetup(X, Y, true_val, total_x_grid_points, num_segments, coord_obs, betas, lambdas)

%start = tic;
[RE] = compute_reflectance("backward", total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y);
%stop = toc(start);
%fprintf("\nIt took %.4f seconds to compute the reflectance.\n\n",stop)

[n1, n2, n3] = size(true_val);
beta = 10^5;
L = -1/2*beta*sqrt(sum(abs(true_val-RE).^2,'all'))/(n1*n2*n3);