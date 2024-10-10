function setup_settings(protein_structure, total_x_grid_points, num_segments, illustrate)

dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
dir = '/zhome/00/b/147112/bioscat/';
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

% Do precomputations
x = linspace(min(X),max(X),total_x_grid_points);
d = x(2)-x(1);
alpha = 2*d;

% Get grid size
n = round(total_x_grid_points/num_segments);
m = num_segments * n;

x = linspace(min(X),max(X),m + 1);
y = interp1(X,Y,x);

segments = cell(num_segments,1);

% Setup segments
for k = 1:num_segments

    segx = x((k-1)*n+1:k*n+1);
    segy = y((k-1)*n+1:k*n+1);
    segments{k} = setup_nanostructures(segx,segy,alpha);

end

save(sprintf('%sData/segments_2D/%s_total_x_grid_points_%d_num_segments_%d.mat',dir,protein_structure,total_x_grid_points,num_segments),"segments")