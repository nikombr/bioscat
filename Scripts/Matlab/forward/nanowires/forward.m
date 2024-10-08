function [C,D,x_int,y_int,x_ext,y_ext] = forward(nanowires)

m = length(nanowires);

% Make vectors ready for shared points
phi = []; x = []; y = []; x_int = []; y_int = []; x_ext = []; y_ext = [];

% Get auxiliary sources and test points
for k = 1:m
    % Particular nanowire
    nw = nanowires{k};
    
    % Combine all points for computation
    phi = [phi; nw.phi];
    x = [x; nw.x];
    y = [y; nw.y];
    x_int = [x_int nw.x_int];
    y_int = [y_int nw.y_int];
    x_ext = [x_ext nw.x_ext];
    y_ext = [y_ext nw.y_ext];
end

coord = struct;
coord_int = struct;
coord_ext = struct;
coord.x = x;
coord.y = y;
coord_int.x = x_int;
coord_int.y = y_int;
coord_ext.x = x_ext;
coord_ext.y = y_ext;

% Compute incident and reflected fields
scenario = 1;
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
b2 = sin(phi) .* (Hx_inc + Hx_ref);
b = [b1; b2];

% Compute scattered fields
Ez_scat = Ez_scat_matrix(coord, coord_int);
Hx_scat = Hx_scat_matrix(coord, coord_int);
Hy_scat = Hy_scat_matrix(coord, coord_int);

% Compute total fields inside the nanowires
Ez_tot_inside = Ez_tot_inside_matrix(coord, coord_ext);
Hx_tot_inside = Hx_tot_inside_matrix(coord, coord_ext);
Hy_tot_inside = Hy_tot_inside_matrix(coord, coord_ext);
    
% Steup matrix
a1 = Ez_scat;
a2 = -Ez_tot_inside;
a3 = -sin(phi) .* Hx_scat       + cos(phi) .* Hy_scat;
a4 =  sin(phi) .* Hx_tot_inside - cos(phi) .* Hy_tot_inside;
A = [a1 a2; a3 a4];


% Solve linear system
c = A\b;

% Save result
C = c(1:length(c)/2);
D = c(length(c)/2+1:end);

% Look at the error of the linear system
error = max(abs(A*c - b));
fprintf("\nThe error from solving the linear system was %.4e\n\n",error);


