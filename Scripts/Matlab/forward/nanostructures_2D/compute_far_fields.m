function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(coord,segments,show_waitbar,scenario)

if nargin < 4
    % Determines the polarisation of the incident and reflected plane wave
    scenario = 1;
end

if nargin < 4
    show_waitbar = true;
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
if show_waitbar
    wb = waitbar(0,'Computing fields...');
    br = 0;
end 
for k = 1:length(segments)
    % Load segment values
    segment    = segments{k};
    coord_int  = segment.coord_int;
    C          = segment.C;

    % Compute scattered fields
    Escat(:,3) = Escat(:,3) + Ez_scat_far_matrix(coord,coord_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_far_matrix(coord,coord_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_far_matrix(coord,coord_int) * C;
    
    if show_waitbar
        br=br+1;
        waitbar(br/(length(segments)),wb);
    end

end
if show_waitbar
    close(wb);
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

% Reshape back to the original shape
Etot  = reshape(Etot,n,m,3);
Htot  = reshape(Htot,n,m,3);
Einc  = reshape(Einc,n,m,3);
Hinc  = reshape(Hinc,n,m,3);
Eref  = reshape(Eref,n,m,3);
Href  = reshape(Href,n,m,3);
Escat = reshape(Escat,n,m,3);
Hscat = reshape(Hscat,n,m,3);

% remove dimensions with 1
if n == 1
    Etot  = squeeze(Etot);
    Htot  = squeeze(Htot);
    Einc  = squeeze(Einc);
    Hinc  = squeeze(Hinc);
    Eref  = squeeze(Eref);
    Href  = squeeze(Href);
    Escat = squeeze(Escat);
    Hscat = squeeze(Hscat);
end
