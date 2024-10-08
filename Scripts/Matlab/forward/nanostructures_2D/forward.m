function segment = forward(segment,scenario,lambda0)

if nargin < 2
    scenario = 1; % If no scenario is chosen, we run the first
end

if nargin < 3
    lambda0 = 325*10^(-9); % Default value of wavelength in free space
end

coord = struct;

% Load vectors of interest
x_top    = segment.x_top;
y_top    = segment.y_top;
n_x      = segment.n_x;
n_y      = segment.n_y;
x_right  = segment.x_right;
y_right  = segment.y_right;
x_bottom = segment.x_bottom;
y_bottom = segment.y_bottom;
x_left   = segment.x_left;
y_left   = segment.y_left;
x_int    = segment.x_int;
y_int    = segment.y_int;
x_ext    = segment.x_ext;
y_ext    = segment.y_ext;

% Store auxiliary points
coord_int = struct;
coord_ext = struct;
coord_int.x = x_int;
coord_int.y = y_int;
coord_ext.x = x_ext;
coord_ext.y = y_ext;

% Top of nanostructure

coord.x = x_top;
coord.y = y_top;

section = 'top';

[b1, b2, a1, a2, a3, a4] = forward_section(coord, coord_int, coord_ext, scenario, lambda0, section, n_x, n_y);

b = [b1; b2];
A = [a1 a2; a3 a4];

% Right side of nanostructure

coord.x = x_right;
coord.y = y_right;

section = 'right';

[b1, b2, a1, a2, a3, a4] = forward_section(coord, coord_int, coord_ext, scenario, lambda0, section);

b = [b; b1; b2];
A = [A; a1 a2; a3 a4];


% Bottom of nanostructure

coord.x = x_bottom;
coord.y = y_bottom;

section = 'bottom';

[b1, b2, a1, a2, a3, a4] = forward_section(coord, coord_int, coord_ext, scenario, lambda0, section);

b = [b; b1; b2];
A = [A; a1 a2; a3 a4];

% Left side of nanostructure

coord.x = x_left;
coord.y = y_left;

section = 'left';

[b1, b2, a1, a2, a3, a4] = forward_section(coord, coord_int, coord_ext, scenario, lambda0, section);

b = [b; b1; b2];
A = [A; a1 a2; a3 a4];

% Solve linear system
c = A\b;

% Save result
segment.C = c(1:length(x_int));
segment.D = c(length(x_int)+1:end);

% Look at the error of the linear system
error = max(abs(A*c - b));
fprintf("\nThe error from solving the linear system was %.4e\n\n",error);


