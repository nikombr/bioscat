%% Compute reflectance and generate synthetic data

clear; close all; clc;
protein_structure    = "Retinin2x2"; % Retinin2x2 demoleus2x2
num_segments = 10; % 1 5 10 20
total_x_grid_points  = 1000;

coord_obs = struct;
r = 3*10^(-2);
phi = linspace(pi/2,2*pi/3,100);
coord_obs.x = cos(phi)*r;
coord_obs.y = sin(phi)*r;
far_field_approximation = false;

lambda0 = 325*10^(-9); % Default value of wavelength in free space

lambdas = linspace(0.5,1.5,200)*lambda0;

betas = linspace(0,pi/2,100);


RE = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas);

% Save clean data (without noise)
filename = sprintf('../../../../Data/reflectance_2D/clean/%s_total_x_grid_points_%d_num_segments_%d.mat',protein_structure,total_x_grid_points,num_segments);
save(filename,'betas','lambdas','RE', 'coord_obs')


%% Add white noise to data

clear; close all; clc;

protein_structure    = "Retinin2x2"; % Retinin2x2 demoleus2x2
num_segments = 1; % 1 5 10 20
total_x_grid_points  = 1000;

for num_segments = [1 5 10 20]
    % Load clean data (without noise)
    filename = sprintf('../../../../Data/reflectance_2D/clean/%s_total_x_grid_points_%d_num_segments_%d.mat',protein_structure,total_x_grid_points,num_segments);
    load(filename)
    
    sigma = 0.01;
    
    RE = RE + randn(size(RE))*sigma;
    RH = RH + randn(size(RH))*sigma;
    
    % Save noisy data
    filename = sprintf('../../../../Data/reflectance_2D/noisy/%s_total_x_grid_points_%d_num_segments_%d.mat',protein_structure,total_x_grid_points,num_segments);
    save(filename,'betas','lambdas','RH','RE', 'coord_obs')

end