function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(coord,coord_int,C,scenario)

if nargin < 4
    % Determines the polarisation of the incident and reflected plane wave
    scenario = 1;
end

% Get coordinates
%X = coord.x;
%Y = coord.y;
x_int = coord_int.x;
y_int = coord_int.y;

% Get size of field to compute
[n,m] = size(coord.x);
coord.x = reshape(coord.x,[],1);
coord.y = reshape(coord.y,[],1);
M = n*m;

% Get the number of seperate computations that have been done
[~, num_computations] = size(C);

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario);
[Eref, Href] = reflected_fields(coord, scenario);

% Compute scattered fields
Escat = zeros(M,3);
Hscat = zeros(M,3);
if num_computations > 1
    coord_temp = struct;
    for k = 1:num_computations
        coord_temp.x = x_int(k,:);
        coord_temp.y = y_int(k,:);
        Escat(:,3) = Escat(:,3) + Ez_scat_far_matrix(coord, coord_temp) * C(:,k);
        Hscat(:,1) = Hscat(:,1) + Hx_scat_far_matrix(coord, coord_temp) * C(:,k);
        Hscat(:,2) = Hscat(:,2) + Hy_scat_far_matrix(coord, coord_temp) * C(:,k);
    end
else
    Escat(:,3) = Escat(:,3) + Ez_scat_far_matrix(coord,coord_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_far_matrix(coord,coord_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_far_matrix(coord,coord_int) * C;
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

% Reshape result
Etot  = reshape(Etot,n,m,3);
Htot  = reshape(Htot,n,m,3);
Einc  = reshape(Einc,n,m,3);
Hinc  = reshape(Hinc,n,m,3);
Eref  = reshape(Eref,n,m,3);
Href  = reshape(Href,n,m,3);
Escat = reshape(Escat,n,m,3);
Hscat = reshape(Hscat,n,m,3);

