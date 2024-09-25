function [Etot, Htot, Escat, Hscat, Einc, Hinc, Eref, Href] = compute_far_fields(X,Y,C,x_int,y_int)

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants();

% Get size of field to compute
[n,m] = size(X);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
M = n*m;

% Get the number of seperate computations that have been done
[~, num_computations] = size(C);

% Compute incident fields
Einc = zeros(M,3);
Einc(:,3) = Ez_inc_vector(X,Y);
Hinc = zeros(M,3);
Hinc(:,1) = Hx_inc_vector(X,Y);

% Compute reflected field
Eref = zeros(M,3);
Eref(:,3) = Ez_ref_vector(X,Y);
Href = zeros(M,3);
Href(:,1) = Hx_ref_vector(X,Y);

% Compute scattered fields
Escat = zeros(M,3);
Hscat = zeros(M,3);
if num_computations > 1
    for k = 1:num_computations
        Escat(:,3) = Escat(:,3) + Ez_scat_far_matrix(X,Y,x_int(k,:),y_int(k,:)) * C(:,k);
        Hscat(:,1) = Hscat(:,1) + Hx_scat_far_matrix(X,Y,x_int(k,:),y_int(k,:)) * C(:,k);
        Hscat(:,2) = Hscat(:,2) + Hy_scat_far_matrix(X,Y,x_int(k,:),y_int(k,:)) * C(:,k);
    end
else
    Escat(:,3) = Escat(:,3) + Ez_scat_far_matrix(X,Y,x_int,y_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_far_matrix(X,Y,x_int,y_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_far_matrix(X,Y,x_int,y_int) * C;
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

% Reshape result
Etot = reshape(Etot,n,m,3);
Htot = reshape(Htot,n,m,3);
Einc = reshape(Einc,n,m,3);
Hinc = reshape(Hinc,n,m,3);
Eref = reshape(Eref,n,m,3);
Href = reshape(Href,n,m,3);
Escat = reshape(Escat,n,m,3);
Hscat = reshape(Hscat,n,m,3);

