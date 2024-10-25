function [X, Y, true_val, args] = nanostructures2DtrueValue(protein_structure, total_x_grid_points, num_segments, data_quality)


filename = sprintf('../../../Data/reflectance_2D/%s/%s_total_x_grid_points_%d_num_segments_%d.mat',data_quality,protein_structure,total_x_grid_points,num_segments);
load(filename,'RE','coord_obs','lambdas','betas')

filename = sprintf("../../../Data/%s_2D.mat",protein_structure);
load(filename,'X','Y')


true_val = RE;

args = {total_x_grid_points, num_segments, coord_obs, betas, lambdas};