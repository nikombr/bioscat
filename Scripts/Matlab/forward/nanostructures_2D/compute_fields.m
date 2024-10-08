function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_fields(coord, segments, far_field_approximation, scenario, lambda0)

if nargin < 4
    % Determines the polarisation of the incident and reflected plane wave
    scenario = 1;
end

if nargin < 5
    lambda0 = 325*10^(-9); % Default value of wavelength in free space
end

% Get size of field to compute
[n,m] = size(coord.x);
coord.x = reshape(coord.x,[],1);
coord.y = reshape(coord.y,[],1);
M = n*m;

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario, lambda0);
[Eref, Href] = reflected_fields(coord, scenario, lambda0);

Escat = zeros(M,3);
Hscat = zeros(M,3);
Etot_inside = zeros(M,3);
Htot_inside = zeros(M,3);
if ~far_field_approximation
    interiors = zeros(M,1);
    num = 2;
end

if canUseGPU
    Escat = gpuArray(Escat);
    Hscat = gpuArray(Hscat);
end

maxtop = -1;

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
    [Escat_matrix, Hscat_matrix] = scattered_fields(coord, coord_int, scenario, lambda0, far_field_approximation);
    if canUseGPU
        tic;
        C = gpuArray(C);
        Escat_matrix = gpuArray(Escat_matrix);
        Hscat_matrix = gpuArray(Hscat_matrix);
        toc;
    end
    if scenario == 1
        Escat(:,3) = Escat(:,3) + Escat_matrix(:,:,3) * C;
        Hscat(:,1) = Hscat(:,1) + Hscat_matrix(:,:,1) * C;
        Hscat(:,2) = Hscat(:,2) + Hscat_matrix(:,:,2) * C;
    elseif scenario ==  2
        Escat(:,1) = Escat(:,1) + Escat_matrix(:,:,1) * C;
        Escat(:,2) = Escat(:,2) + Escat_matrix(:,:,2) * C;
        Hscat(:,3) = Hscat(:,3) + Hscat_matrix(:,:,3) * C;
    end
    %Escat(:,3) = Escat(:,3) + Ez_scat_matrix(coord, coord_int) * C;
    %Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(coord, coord_int) * C;
    %Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(coord, coord_int) * C;

    br=br+1;
    waitbar(br/(num*length(segments)),wb);
    
    hej = max(segment.y_top) < min(coord.y)
    maxtop = max(maxtop,max(segment.y_top));
    if ~far_field_approximation && max(segment.y_top) < min(coord.y)
        % Find the interior of the nanostructures
        f = interp1(segment.x_top,segment.y_top,coord.x);
        interior = coord.y < f;
        interiors = interiors + interior;
        
        % Get fieds inside the nanostructure
        E_inside = zeros(M,3);
        H_inside = zeros(M,3);

        [Etot_inside_matrix, Htot_inside_matrix] = interior_fields(coord, coord_ext, scenario, lambda0);

        if scenario == 1
            E_inside(:,3) = Etot_inside_matrix(:,:,3) * D;
            H_inside(:,1) = Htot_inside_matrix(:,:,1) * D;
            H_inside(:,2) = Htot_inside_matrix(:,:,2) * D;
        elseif scenario ==  2
            E_inside(:,1) = Etot_inside_matrix(:,:,1) * D;
            E_inside(:,2) = Etot_inside_matrix(:,:,2) * D;
            H_inside(:,3) = Htot_inside_matrix(:,:,3) * D;
        end
        
        %E_inside(:,3) = Ez_tot_inside_matrix(coord, coord_ext) * D;
        %H_inside(:,1) = Hx_tot_inside_matrix(coord, coord_ext) * D;
        %H_inside(:,2) = Hy_tot_inside_matrix(coord, coord_ext) * D;
        
        % Remove field outside of nanostructures
        E_inside(~interior,:) = 0;
        H_inside(~interior,:) = 0;
        
        % Add total field inside nanostructures
        Etot_inside = Etot_inside + E_inside;
        Htot_inside = Htot_inside + H_inside;
    end

    br=br+1;
    waitbar(br/(num*length(segments)),wb);
end
close(wb);

if canUseGPU
    Escat = gather(Escat);
    Hscat = gather(Hscat);
end

if maxtop < min(coord.y)
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
else
    % Get the total fields
    Etot = Einc + Eref + Escat;
    Htot = Hinc + Href + Hscat;
end



% Reshape back to the original shape
Etot = reshape(Etot,n,m,3);
Htot = reshape(Htot,n,m,3);
Einc = reshape(Einc,n,m,3);
Hinc = reshape(Hinc,n,m,3);
Eref = reshape(Eref,n,m,3);
Href = reshape(Href,n,m,3);
Escat = reshape(Escat,n,m,3);
Hscat = reshape(Hscat,n,m,3);

