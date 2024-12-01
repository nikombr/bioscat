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
total_grid_points = 300;
obs_grid = 500;
beta = 90;
protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2
filename = sprintf("%s/Results/forward/%s/near/fields_beta_%d_lambda_%d_num_segments_1_total_grid_points_%d_obs_grid_%d.mat",dir,protein_structure,beta,lambda,total_grid_points,obs_grid);
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

% Set nanostructure (have done this manually in comsol as we only have two nanostructures)
%surface = load(sprintf('%s/Data/comsol/surface_%s.txt',dir,protein_structure))
%surface_interpolation = model.func('int1')
%surface_interpolation.set('table', cellfun( @(x) num2str(x, '%0.10g'), num2cell( surface ), 'UniformOutput', false));

% Set the parameter in the COMSOL model
param_name = 'lambda_0';    % Name of the parameter in the COMSOL model
param_value = lambda*1e-9;          % New value for the parameter
model.param.set(param_name, param_value);
param_name = 'beta';    % Name of the parameter in the COMSOL model
param_value = beta*pi/180;          % New value for the parameter
model.param.set(param_name, param_value);

%emw = model.physics('emw');

% Access the feature that defines the electric field
%formulation = emw.prop('BackgroundField') % Replace 'userDefinedField' with the actual feature tag
%disp(emw.prop.tags());
% Update the field components
%disp(formulation.tags());
%formulation.get()
%formulation.set('field', '{exp(1j*k_0*y), 0, 0}');
%formulation.set('Ey', '0');
%formulation.set('Ez', '0');

% Solve
model.sol('sol1').runAll;

vars = {'x', 'y', 'z'};

E_comsol = E_tot*0;
H_comsol = E_tot*0;

% Get solution
for i = 1:3
    result = mpheval(model,sprintf('emw.E%s',vars{i}));
    values = result.d1;  % Field values
    coord = result.p;  % Coordonates
    x_mesh = coord(1, :);  % x-coordinates of the mesh
    y_mesh = coord(2, :);  % y-coordinates of the mesh
    E_comsol(:,:,i) = griddata(x_mesh, y_mesh, values, Xmesh, Ymesh, 'linear'); % Interpolate the field values at the desired points
end

for i = 1:3
    result = mpheval(model,sprintf('emw.H%s',vars{i}));
    values = result.d1;  % Field values
    coord = result.p;  % Coordonates
    x_mesh = coord(1, :);  % x-coordinates of the mesh
    y_mesh = coord(2, :);  % y-coordinates of the mesh
    H_comsol(:,:,i) = griddata(x_mesh, y_mesh, values, Xmesh, Ymesh, 'linear'); % Interpolate the field values at the desired points
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

