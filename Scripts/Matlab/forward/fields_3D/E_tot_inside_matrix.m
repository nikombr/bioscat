function E_tot_inside = E_tot_inside_matrix(x,y,z,x_ext,y_ext,z_ext)

% Allocate vectors
E_tot_inside = zeros(length(x),length(x_ext),3);

% Translation of coordinate system to have origo at (x_ext,y_ext,z_ext)
X = x - x_ext;
Y = y - y_ext;
Z = z - z_ext;

% Get fields
[EHx, EHy, EHz] = EH(X, Y, Z);

% Compute fields in the Cartesian coordinate system
E_tot_inside(:,:,1) = EHx;
E_tot_inside(:,:,2) = EHy;
E_tot_inside(:,:,3) = EHz;