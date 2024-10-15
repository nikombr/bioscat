function setup_settings(protein_structure, total_x_grid_points, num_segments, illustrate)

dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
%dir = '/zhome/00/b/147112/bioscat/';
if nargin < 4
    illustrate = false;
end

load(sprintf("%sData/%s_2D.mat",dir,protein_structure))

if illustrate
    figure('Renderer', 'painters', 'Position', [400 400 800 250]);
    plot(X,Y,'LineWidth',1.5)
    xlim([min(X),max(X)])
    grid on
    hold on
end

segments = setup_segments(X,Y,num_segments,total_x_grid_points);

save(sprintf('%sData/segments_2D/%s_total_x_grid_points_%d_num_segments_%d.mat',dir,protein_structure,total_x_grid_points,num_segments),"segments")