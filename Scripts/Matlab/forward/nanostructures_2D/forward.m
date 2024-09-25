function segment = forward_nanostructures_2D(segment)

addpath('compute_field_vectors')

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, k1] = load_constants_nanostructures_2D();

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

%% Top of nanostructure

x = x_top;
y = y_top;

% Compute incident fields
Ez_inc = Ez_inc_vector(x,y);
Hx_inc = Hx_inc_vector(x,y);

% Compute reflected fields
Ez_ref = Ez_ref_vector(x,y);
Hx_ref = Hx_ref_vector(x,y);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = n_y .* (Hx_inc + Hx_ref);  % z component
b = [b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(x,y,x_int,y_int);
Hx_scat = Hx_scat_matrix(x,y,x_int,y_int);
Hy_scat = Hy_scat_matrix(x,y,x_int,y_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(x,y,x_ext,y_ext);
Hx_tot_inside = Hx_tot_inside_matrix(x,y,x_ext,y_ext);
Hy_tot_inside = Hy_tot_inside_matrix(x,y,x_ext,y_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 =   n_x .* Hy_scat - n_y .* Hx_scat;
a4 = - n_x .* Hy_tot_inside + n_y .* Hx_tot_inside;
A = [a1 a2; a3 a4];

%% Right side of nanostructure

x = x_right;
y = y_right;

% Compute incident fields
Ez_inc = Ez_inc_vector(x,y);

% Compute reflected fields
Ez_ref = Ez_ref_vector(x,y);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = 0*Ez_inc;
b = [b; b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(x,y,x_int,y_int);
Hy_scat = Hy_scat_matrix(x,y,x_int,y_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(x,y,x_ext,y_ext);
Hy_tot_inside = Hy_tot_inside_matrix(x,y,x_ext,y_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 = Hy_scat;
a4 = -Hy_tot_inside;
A = [A; a1 a2; a3 a4];


%% Bottom of nanostructure

x = x_bottom;
y = y_bottom;

% Compute incident fields
Ez_inc = Ez_inc_vector(x,y);
Hx_inc = Hx_inc_vector(x,y);

% Compute reflected fields
Ez_ref = Ez_ref_vector(x,y);
Hx_ref = Hx_ref_vector(x,y);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = Hx_inc + Hx_ref;
b = [b; b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(x,y,x_int,y_int);
Hx_scat = Hx_scat_matrix(x,y,x_int,y_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(x,y,x_ext,y_ext);
Hx_tot_inside = Hx_tot_inside_matrix(x,y,x_ext,y_ext);
    
% Setup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 = Hx_scat;
a4 = -Hx_tot_inside;
A = [A; a1 a2; a3 a4];

%% Left side of nanostructure

x = x_left;
y = y_left;

% Compute incident fields
Ez_inc = Ez_inc_vector(x,y);

% Compute reflected fields
Ez_ref = Ez_ref_vector(x,y);
    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = 0*Ez_inc;
b = [b; b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(x,y,x_int,y_int);
Hy_scat = Hy_scat_matrix(x,y,x_int,y_int);

% Compute total fields inside the nanostructure
Ez_tot_inside = Ez_tot_inside_matrix(x,y,x_ext,y_ext);
Hy_tot_inside = Hy_tot_inside_matrix(x,y,x_ext,y_ext);
    
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


