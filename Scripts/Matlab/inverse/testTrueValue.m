function [X, Y, true_val, args] = testTrueValue(protein_structure, total_x_grid_points, num_segments, data_quality)


filename = sprintf("../../../Data/%s_2D.mat",protein_structure);
load(filename,'X','Y')


true_val = Y;
args = {};