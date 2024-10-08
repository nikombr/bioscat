function H_tot_inside = H_tot_inside_matrix(x,y,z,x_ext,y_ext,z_ext)

% Allocate vectors
H_tot_inside = zeros(length(x),length(x_ext),3);

% Translation of coordinate system to have origo at (x_ext,y_ext,z_ext)
X = x - x_ext;
Y = y - y_ext;
Z = z - z_ext;

% Compute fields
[HHx, HHy, HHz] = HH(X,Y,Z);

% Compute fields in the Cartesian coordinate system
H_tot_inside(:,:,1) = HHx;
H_tot_inside(:,:,2) = HHy;
H_tot_inside(:,:,3) = HHz;