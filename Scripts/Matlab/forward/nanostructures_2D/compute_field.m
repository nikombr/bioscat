function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_2D_field_nanostructures_2D(X,Y,segments,C,D,x_int,y_int,x_ext,y_ext)

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, k1] = load_constants_nanostructures_2D();

% Get size of field to compute
[n,m] = size(X);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
M = n*m;

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

Escat = zeros(M,3);
Hscat = zeros(M,3);
Etot_inside = zeros(M,3);
Htot_inside = zeros(M,3);
interiors = zeros(M,1);

wb = waitbar(0,'Computing fields...');
br = 0;
for k = 1:length(segments)
    % Load segment values
    segment = segments{k};
    x_int  = segment.x_int;
    y_int  = segment.y_int;
    x_ext  = segment.x_ext;
    y_ext  = segment.y_ext;
    C      = segment.C;
    D      = segment.D;

    % Compute scattered fields
    Escat(:,3) = Escat(:,3) + Ez_scat_matrix(X,Y,x_int,y_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(X,Y,x_int,y_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(X,Y,x_int,y_int) * C;

    br=br+1;
    waitbar(br/(2*length(segments)),wb);
    
    % Find the interior of the nanostructures
    f = interp1(segment.x_top,segment.y_top,X);
    interior = Y < f;
    interiors = interiors + interior;
    
    % Get fieds inside the nanostructure
    E_inside = zeros(M,3);
    H_inside = zeros(M,3);
    
    E_inside(:,3) = Ez_tot_inside_matrix(X,Y,x_ext,y_ext) * D;
    H_inside(:,1) = Hx_tot_inside_matrix(X,Y,x_ext,y_ext) * D;
    H_inside(:,2) = Hy_tot_inside_matrix(X,Y,x_ext,y_ext) * D;
    
    % Remove field outside of nanostructures
    E_inside(~interior,:) = 0;
    H_inside(~interior,:) = 0;
    
    % Add total field inside nanostructures
    Etot_inside = Etot_inside + E_inside;
    Htot_inside = Htot_inside + H_inside;

    br=br+1;
    waitbar(br/(2*length(segments)),wb);
end
close(wb);

% Remove the fields in the interior of nanostructures
remove = interiors > 0;
Escat(remove, :) = 0;
Eref(remove,  :) = 0;
Einc(remove,  :) = 0;
Hscat(remove, :) = 0;
Href(remove,  :) = 0;
Hinc(remove,  :) = 0;

% Get the total fields
Etot = Etot_inside + Einc + Eref + Escat;
Htot = Htot_inside + Hinc + Href + Hscat;

% Reshape back to the original shape
Etot = reshape(Etot,n,m,3);
Htot = reshape(Htot,n,m,3);
Einc = reshape(Einc,n,m,3);
Hinc = reshape(Hinc,n,m,3);
Eref = reshape(Eref,n,m,3);
Href = reshape(Href,n,m,3);
Escat = reshape(Escat,n,m,3);
Hscat = reshape(Hscat,n,m,3);







