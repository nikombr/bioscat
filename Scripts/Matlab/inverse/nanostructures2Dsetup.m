function L = nanostructures2Dsetup(X, Y, true_val, total_x_grid_points, num_segments, coord_obs, betas, lambdas)

start = tic;
[RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y);
stop = toc(start);
fprintf("\nIt took %.4f seconds to compute the reflectance.\n\n",stop)

L = -sqrt(sum(abs(true_val-RE).^2,'all'))/2;