
clear; close all; clc;

protein_structure = 'demoleus2x2'; % demoleus2x2 Retinin2x2

scenario = 1;

if scenario == 1
    beta = 0;
elseif scenario == 2
    beta = 90;
end
lambda=325; % nm
total_grid_points = 300;

filename = sprintf("../Results/forward/%s/far_field_pattern/fields_beta_%d_lambda_%d_num_segments_1_total_grid_points_%d.mat",protein_structure,beta,lambda,total_grid_points);
load(filename)

% Define norm
normfunc = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

% Compute far-field pattern
PE = normfunc(abs(E_scat+E_ref));
PH = normfunc(abs(H_scat));
phi = linspace(0,pi,200);

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(1,1,'TileSpacing','compact');

nexttile
polarplot(phi,PE,'-','linewidth',1.5)
thetalim([0 180])
if strcmp(protein_structure,'Retinin2x2')
    %rlim([0,0.025])
elseif strcmp(protein_structure,'demoleus2x2')
    %rlim([0,0.03])
end
thetaticks([0 30 60 90 120 150 180])
thetaticklabels({'$0$','$\pi/6$','$\pi/3$','$\pi/2$','$2\pi/3$','$5\pi/6$','$\pi$'})
%text(3*pi/2,0.005,'Intensity')

l = legend('numcolumns',2,'FontSize',14,'box','off');
l.Layout.Tile = 'south';

destination = sprintf('../Illustrations/nanostructures/%s_far_field_patterns.png',protein_structure);
exportgraphics(gcf,destination,'Resolution',300);