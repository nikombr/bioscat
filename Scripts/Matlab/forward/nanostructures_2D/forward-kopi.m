function segment = forward(segment,scenario)

if nargin < 2
    scenario = 1;
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

coord_int = struct;
coord_ext = struct;
coord_int.x = x_int;
coord_int.y = y_int;
coord_ext.x = x_ext;
coord_ext.y = y_ext;

%% Top of nanostructure

coord.x = x_top;
coord.y = y_top;

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario);
[Eref, Href] = reflected_fields(coord, scenario);

% Compute incident fields
Ez_inc = Einc(:,3);
Hx_inc = Hinc(:,1);

% Compute reflected fields
Ez_ref = Eref(:,3);
Hx_ref = Href(:,1);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = n_y .* (Hx_inc + Hx_ref);  % z component
b = [b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(coord, coord_int);
Hx_scat = Hx_scat_matrix(coord, coord_int);
Hy_scat = Hy_scat_matrix(coord, coord_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(coord, coord_ext);
Hx_tot_inside = Hx_tot_inside_matrix(coord, coord_ext);
Hy_tot_inside = Hy_tot_inside_matrix(coord, coord_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 =   n_x .* Hy_scat - n_y .* Hx_scat;
a4 = - n_x .* Hy_tot_inside + n_y .* Hx_tot_inside;
A = [a1 a2; a3 a4];

%% Right side of nanostructure

coord.x = x_right;
coord.y = y_right;

% Compute incident and reflected fields
[Einc, ~] =  incident_fields(coord, scenario);
[Eref, ~] = reflected_fields(coord, scenario);

% Compute incident fields
Ez_inc = Einc(:,3);

% Compute reflected fields
Ez_ref = Eref(:,3);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = 0*Ez_inc;
b = [b; b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(coord, coord_int);
Hy_scat = Hy_scat_matrix(coord, coord_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(coord, coord_ext);
Hy_tot_inside = Hy_tot_inside_matrix(coord, coord_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 = Hy_scat;
a4 = -Hy_tot_inside;
A = [A; a1 a2; a3 a4];


%% Bottom of nanostructure

coord.x = x_bottom;
coord.y = y_bottom;

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario);
[Eref, Href] = reflected_fields(coord, scenario);

% Compute incident fields
Ez_inc = Einc(:,3);
Hx_inc = Hinc(:,1);

% Compute reflected fields
Ez_ref = Eref(:,3);
Hx_ref = Href(:,1);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = Hx_inc + Hx_ref;
b = [b; b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(coord, coord_int);
Hx_scat = Hx_scat_matrix(coord, coord_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(coord, coord_ext);
Hx_tot_inside = Hx_tot_inside_matrix(coord, coord_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 = Hx_scat;
a4 = -Hx_tot_inside;
A = [A; a1 a2; a3 a4];

%% Left side of nanostructure

coord.x = x_left;
coord.y = y_left;

% Compute incident and reflected fields
[Einc, ~] =  incident_fields(coord, scenario);
[Eref, ~] = reflected_fields(coord, scenario);

% Compute incident fields
Ez_inc = Einc(:,3);

% Compute reflected fields
Ez_ref = Eref(:,3);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = 0*Ez_inc;
b = [b; b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(coord, coord_int);
Hy_scat = Hy_scat_matrix(coord, coord_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(coord, coord_ext);
Hy_tot_inside = Hy_tot_inside_matrix(coord, coord_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 = Hy_scat;
a4 = -Hy_tot_inside;
A = [A; a1 a2; a3 a4];
%% Finalize

% Solve linear system
c = A\b;

% Save result
segment.C = c(1:length(x_int));
segment.D = c(length(x_int)+1:end);

% Look at the error of the linear system
error = max(abs(A*c - b));
fprintf("\nThe error from solving the linear system was %.4e\n\n",error);


