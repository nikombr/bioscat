function E_scat = E_scat_matrix(x,y,z,x_int,y_int,z_int)
% Allocate vectors
E_scat = zeros(length(x),length(x_int),3);

% Translation of coordinate system to have origo at (x_int,y_int,z_int)
X     = x - x_int;
Y     = y - y_int;
Y_ref = y + y_int;
Z     = z - z_int;

% Get fields
[EHx,     EHy,     EHz]     = EH(X, Y,     Z);
[EHx_ref, EHy_ref, EHz_ref] = EH(X, Y_ref, Z);

% Compute fields in the Cartesian coordinate system
E_scat(:,:,1) = EHx + Gamma_r * EHx_ref;
E_scat(:,:,2) = EHy + Gamma_r * EHy_ref;
E_scat(:,:,3) = EHz + Gamma_r * EHz_ref;