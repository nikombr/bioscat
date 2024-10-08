function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_fields(coord,coord_int,coord_ext,nanowires,C,D,scenario)

if nargin < 7
    % Determines the polarisation of the incident and reflected plane wave
    scenario = 1;
end

x_int = coord_int.x;
y_int = coord_int.y;
x_ext = coord_ext.x;
y_ext = coord_ext.y;

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
        Escat(:,3) = Escat(:,3) + Ez_scat_matrix(coord, coord_temp) * C(:,k);
        Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(coord, coord_temp) * C(:,k);
        Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(coord, coord_temp) * C(:,k);
    end
else
    Escat(:,3) = Escat(:,3) + Ez_scat_matrix(coord,coord_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(coord,coord_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(coord,coord_int) * C;
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

% Find the interior of the nanowires
numerical = @(x,y) sqrt(x.^2+y.^2);
find_interior = zeros(m*n,1);
for j = 1:length(nanowires)
    nw = nanowires{j};
    dist = numerical(coord.x - nw.xc, coord.y - nw.r);
    find_interior = find_interior + (dist < nw.r);
end

% Remove the fields form the interior of the nanowires
Etot(find_interior > 0,:) = 0;
Htot(find_interior > 0,:) = 0;


% Compute total fields inside nanowires
if num_computations > 1
    coord_temp = struct;
    for k = 1:num_computations
        Etot_inside = zeros(M,3);
        Htot_inside = zeros(M,3);
        coord_temp.x = x_ext(k,:);
        coord_temp.y = y_ext(k,:);
        Etot_inside(:,3) = Ez_tot_inside_matrix(coord, coord_temp) * D(:,k);
        Htot_inside(:,1) = Hx_tot_inside_matrix(coord, coord_temp) * D(:,k);
        Htot_inside(:,2) = Hy_tot_inside_matrix(coord, coord_temp) * D(:,k);
    
        % Remove field outside of nanowires
        Etot_inside(~find_interior,:) = 0;
        Htot_inside(~find_interior,:) = 0;
        
        % Add total field inside nanowires
        Etot = Etot + Etot_inside;
        Htot = Htot + Htot_inside;
    end
else
    Etot_inside = zeros(M,3);
    Htot_inside = zeros(M,3);

    Etot_inside(:,3) = Ez_tot_inside_matrix(coord, coord_ext) * D;
    Htot_inside(:,1) = Hx_tot_inside_matrix(coord, coord_ext) * D;
    Htot_inside(:,2) = Hy_tot_inside_matrix(coord, coord_ext) * D;

    % Remove field outside of nanowires
    Etot_inside(~find_interior,:) = 0;
    Htot_inside(~find_interior,:) = 0;
    
    % Add total field inside nanowires
    Etot = Etot + Etot_inside;
    Htot = Htot + Htot_inside;

end

Etot  = reshape(Etot,n,m,3);
Htot  = reshape(Htot,n,m,3);
Einc  = reshape(Einc,n,m,3);
Hinc  = reshape(Hinc,n,m,3);
Eref  = reshape(Eref,n,m,3);
Href  = reshape(Href,n,m,3);
Escat = reshape(Escat,n,m,3);
Hscat = reshape(Hscat,n,m,3);

