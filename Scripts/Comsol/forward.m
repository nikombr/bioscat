close all force
clear all
clc
beep off
format compact
format long

% Directory data and results
dir = "../../../../../../../work3/s194146/bioscatdata";

% Load my solution
lambda=325; % nm
total_grid_points = 1000;
obs_grid = 300;
protein_structure = 'Retinin2x2'; % Retinin2x2 demoleus2x2
for beta = [0 30 60 90]
    beta
filename = sprintf("%s/Results/forward/comsol_test/%s/near/fields_beta_%d_lambda_%d_num_segments_1_total_grid_points_%d_obs_grid_%d.mat",dir,protein_structure,beta,lambda,total_grid_points,obs_grid);
load(filename)
E_tot = E_scat + E_inc + E_int;
H_tot = H_scat + H_inc + H_int;

% Get mesh for points
[Xmesh, Ymesh] = meshgrid(x,y);
coords = [Xmesh(:), Ymesh(:)]';     % Convert to a list of points

% Load modules
import com.comsol.model.*
import com.comsol.model.util.*

% Load Comsol model
model = mphload(sprintf('%s/../forward61_%s.mph',dir,protein_structure));

% Set the parameter in the COMSOL model
param_name = 'lambda_0';    % Name of the parameter in the COMSOL model
param_value = lambda*1e-9;          % New value for the parameter
model.param.set(param_name, param_value);
param_name = 'beta';    % Name of the parameter in the COMSOL model
param_value = beta*pi/180;          % New value for the parameter
model.param.set(param_name, param_value);

% Solve
model.sol('sol1').runAll;

vars = {'x', 'y', 'z'};

E_comsol = E_tot*0;
H_comsol = E_tot*0;

for i = 1:3
    result = mphinterp(model, sprintf('emw.E%s',vars{i}), 'coord', coords);
    E_comsol(:,:,i) = reshape(result, [obs_grid obs_grid]);
end

for i = 1:3
    result = mphinterp(model, sprintf('emw.H%s',vars{i}), 'coord', coords);
    H_comsol(:,:,i) = reshape(result, [obs_grid obs_grid]);
end

k0 = 2 * pi / param_value;  % Wave number (lambda is the wavelength)
E_background = exp(1j * k0 * Ymesh);  % Background field

E_diff = E_comsol - E_tot;
H_diff = H_comsol - H_tot;

% Plot the interpolated field
figure('Renderer', 'painters', 'Position', [400 400 1500 1500]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:3
    nexttile;
    imagesc(x, y, abs(E_comsol(:,:,i)));  % Absolute value of Ez for visualization
    title('Comsol');
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel(sprintf('|E%s|',vars{i}));
    axis square
    colorbar
    %clim([0,2])
    ax = gca;              % Get the current axes
    ax.CLim = [0 2.5]; % Set the color limits
    set(gca,'YDir','normal')

    nexttile;
    imagesc(x, y, abs(E_tot(:,:,i)));  % Absolute value of Ez for visualization
    axis square
    colorbar
    %clim([0,2])
    ax = gca;              % Get the current axes
    ax.CLim = [0 2.5]; % Set the color limits
    set(gca,'YDir','normal')

    nexttile;
    imagesc(x, y, abs(E_diff(:,:,i))./abs(E_comsol(:,:,i)+1e-20)); 
    axis square
    colorbar
    ax = gca;              % Get the current axes
    %ax.CLim = [0 0.5]; % Set the color limits
    set(gca,'YDir','normal')
    set(gca,'colorscale','log')
end

exportgraphics(gcf,sprintf('plot_%s_E_beta_%d.png',protein_structure,beta),'Resolution',300);

% Plot the interpolated field
figure('Renderer', 'painters', 'Position', [400 400 1500 1500]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:3
    nexttile;
    imagesc(x, y, abs(H_comsol(:,:,i)));  % Absolute value of Ez for visualization
    title('Comsol');
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel(sprintf('|E%s|',vars{i}));
    axis square
    colorbar
    %clim([0,2])
    ax = gca;              % Get the current axes
    %ax.CLim = [0 2.5]; % Set the color limits
    set(gca,'YDir','normal')

    nexttile;
    imagesc(x, y, abs(H_tot(:,:,i)));  % Absolute value of Ez for visualization
    axis square
    colorbar
    %clim([0,2])
    ax = gca;              % Get the current axes
    %ax.CLim = [0 2.5]; % Set the color limits
    set(gca,'YDir','normal')

    nexttile;
    imagesc(x, y, abs(H_diff(:,:,i))./abs(H_comsol(:,:,i)+1e-20)); 
    axis square
    colorbar
    ax = gca;              % Get the current axes
    %ax.CLim = [0 0.5]; % Set the color limits
    set(gca,'YDir','normal')
    set(gca,'colorscale','log')
end

exportgraphics(gcf,sprintf('plot_%s_H_beta_%d.png',protein_structure,beta),'Resolution',300);

filename = sprintf("%s/Data/comsol/%s/near/fields_beta_%d_lambda_%d_total_grid_points_%d_obs_grid_%d.mat",dir,protein_structure,beta,lambda,total_grid_points,obs_grid);
save(filename, 'E_comsol', 'H_comsol', 'x', 'y')

end