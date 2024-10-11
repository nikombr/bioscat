function segments = setup_segments(X,Y,num_segments,total_x_grid_points)


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