function H_scat = H_scat_matrix(x,y,z,x_int,y_int,z_int)

% Allocate
H_scat = zeros(length(x),length(x_int),3);

% Translation of coordinate system to have origo at (x_int,y_int,z_int)
X     = x - x_int;
Y     = y - y_int;
Y_ref = y + y_int;
Z     = z - z_int;

% Get fields
[HHx,     HHy,     HHz]     = HH(X, Y,     Z);
[HHx_ref, HHy_ref, HHz_ref] = HH(X, Y_ref, Z);

% Compute fields in the Cartesian coordinate system
H_scat(:,:,1) = HHx + Gamma_r * HHx_ref;
H_scat(:,:,2) = HHy + Gamma_r * HHy_ref;
H_scat(:,:,3) = HHz + Gamma_r * HHz_ref;