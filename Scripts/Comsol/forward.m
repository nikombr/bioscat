close all force
clear all
clc
beep off
format compact
format long

% Load my solution
lambda=325; % nm
total_grid_points = 300;
beta = 0;
protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2
filename = sprintf("../../Results/forward/2D/%s/fields_%d_lambda_%d_num_segments_1_total_grid_points_%d.mat",protein_structure,beta,lambda,total_grid_points);
load(filename)

% Get mesh for points
[Xmesh, Ymesh] = meshgrid(x,y);
coords = [Xmesh(:), Ymesh(:)]';     % Convert to a list of points

% Load nanostructure
xs = load("../../Data/nanostructures/demoleus2x2_2D_x_300.txt");
fs = load("../../Data/nanostructures/demoleus2x2_2D_f_300.txt");

% Load modules
import com.comsol.model.*
import com.comsol.model.util.*

figure, plot(xs,fs), grid on, axis equal
surface = [xs' fs'];
save('surface.txt','surface','-ascii')

% Load Comsol model
model = mphload('forward61.mph');

% Set the parameter in the COMSOL model
param_name = 'lambda_0';    % Name of the parameter in the COMSOL model
param_value = lambda*1e-9;          % New value for the parameter
model.param.set(param_name, param_value);

% Solve
model.sol('sol1').runAll;

% Get solution
result = mpheval(model,'emw.Ez')
Ez_values = result.d1;  % Field values
coord = result.p;  % Coordonates
x_mesh = coord(1, :);  % x-coordinates of the mesh
y_mesh = coord(2, :);  % y-coordinates of the mesh
% Interpolate the field values at the desired points
Ez_interp = griddata(x_mesh, y_mesh, Ez_values, Xmesh, Ymesh, 'linear');

k0 = 2 * pi / param_value;  % Wave number (lambda is the wavelength)
E_background = exp(1j * k0 * Ymesh);  % Background field

% Plot the interpolated field
figure('Renderer', 'painters', 'Position', [400 400 1500 550]);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
imagesc(x, y, abs(Ez_interp-E_background));  % Absolute value of Ez for visualization
title('Interpolated emw.Ez');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('|Ez|');
grid on;
axis square
colorbar
%clim([0,2])
ax = gca;              % Get the current axes
ax.CLim = [0 1.8]; % Set the color limits

nexttile;
Ez_tot = E_scat(:,:,3)+E_inc(:,:,3)+E_ref(:,:,3);
%Ez_tot = E_scat(:,:,3)+E_inc(:,:,3);
imagesc(x, y, abs(Ez_tot-E_inc(:,:,3)));  % Absolute value of Ez for visualization
axis square
colorbar
%clim([0,2])
ax = gca;              % Get the current axes
ax.CLim = [0 1.8]; % Set the color limits

nexttile;
imagesc(x, y, abs(Ez_tot-Ez_interp)); 
axis square
colorbar

exportgraphics(gcf,'hmm.png','Resolution',300);