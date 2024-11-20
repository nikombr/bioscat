close all force
clear all
clc
beep off
format compact
format long

% Load my solution
lambda=325; % nm
total_grid_points = 1000;
beta = 0;
protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2
filename = sprintf("../../Results/forward/%s/near/fields_beta_%d_lambda_%d_num_segments_1_total_grid_points_%d.mat",protein_structure,beta,lambda,total_grid_points);
load(filename)

% Get mesh for points
[Xmesh, Ymesh] = meshgrid(x,y);
coords = [Xmesh(:), Ymesh(:)]';     % Convert to a list of points

% Load nanostructure
xs = load("../../Data/nanostructures/2D/demoleus2x2_x_300.txt");
fs = load("../../Data/nanostructures/2D/demoleus2x2_f_300.txt");

% Load modules
import com.comsol.model.*
import com.comsol.model.util.*

% Set location for temporary files
-tmpdir $HOME/work3/s194146/comsoltmp
-recoverydir $HOME/work3/s194146/comsolrecovery

figure, plot(xs,fs), grid on, axis equal
surface = [xs' fs'];
save('surface.txt','surface','-ascii')

% Load Comsol model
model = mphload('forward61_2.mph');

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
Ez_comsol = Ez_interp-E_background;
<<<<<<< HEAD
imagesc(x, y, abs(Ez_comsol));  % Absolute value of Ez for visualization
=======
imagesc(x, y, abs(Ez_comsol).^2);  % Absolute value of Ez for visualization
>>>>>>> 6ef0007e1ccbb9c08ab322f0c9352dee59aea4cf
title('Interpolated emw.Ez');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('|Ez|');
axis square
colorbar
%clim([0,2])
ax = gca;              % Get the current axes
<<<<<<< HEAD
ax.CLim = [0 0.5]; % Set the color limits

nexttile;
Ez_me = E_scat(:,:,3) + E_ref(:,:,3);
%Ez_tot = E_scat(:,:,3)+E_inc(:,:,3);
imagesc(x, y, abs(Ez_me));  % Absolute value of Ez for visualization
=======
ax.CLim = [0 0.4]; % Set the color limits

nexttile;
Ez_me = E_scat(:,:,3)+E_ref(:,:,3);
%Ez_tot = E_scat(:,:,3)+E_inc(:,:,3);
imagesc(x, y, abs(Ez_me).^2);  % Absolute value of Ez for visualization
>>>>>>> 6ef0007e1ccbb9c08ab322f0c9352dee59aea4cf
axis square
colorbar
%clim([0,2])
ax = gca;              % Get the current axes
<<<<<<< HEAD
ax.CLim = [0 0.5]; % Set the color limits

nexttile;
imagesc(x, y, abs(Ez_comsol-Ez_me)); 
=======
ax.CLim = [0 0.4]; % Set the color limits

nexttile;
imagesc(x, y, abs(Ez_comsol-Ez_me).^2); 
>>>>>>> 6ef0007e1ccbb9c08ab322f0c9352dee59aea4cf
axis square
colorbar
ax = gca;              % Get the current axes
ax.CLim = [0 0.5]; % Set the color limits

exportgraphics(gcf,'hmm.png','Resolution',300);