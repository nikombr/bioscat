function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_fields(coord,segments,scenario)

if nargin < 4
    % Determines the polarisation of the incident and reflected plane wave
    scenario = 1;
end


% Get size of field to compute
[n,m] = size(coord.x);
coord.x = reshape(coord.x,[],1);
coord.y = reshape(coord.y,[],1);
M = n*m;

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario);
[Eref, Href] = reflected_fields(coord, scenario);

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
    coord_int  = segment.coord_int;
    coord_ext  = segment.coord_ext;
    C      = segment.C;
    D      = segment.D;

    % Compute scattered fields
    Escat(:,3) = Escat(:,3) + Ez_scat_matrix(coord, coord_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(coord, coord_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(coord, coord_int) * C;

    br=br+1;
    waitbar(br/(2*length(segments)),wb);
    
    % Find the interior of the nanostructures
    f = interp1(segment.x_top,segment.y_top,coord.x);
    interior = coord.y < f;
    interiors = interiors + interior;
    
    % Get fieds inside the nanostructure
    E_inside = zeros(M,3);
    H_inside = zeros(M,3);
    
    E_inside(:,3) = Ez_tot_inside_matrix(coord, coord_ext) * D;
    H_inside(:,1) = Hx_tot_inside_matrix(coord, coord_ext) * D;
    H_inside(:,2) = Hy_tot_inside_matrix(coord, coord_ext) * D;
    
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

