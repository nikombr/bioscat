function [E_error, H_error] = transmission_error(coord_int,coord_ext,C,D,nanowires,n)

scenario = 1;

% Get the number of seperate computations that have been done
[~, num_computations] = size(C);

% n is number of points uniformly distributed over the circles, if it is
% not provided, the error is computed in the test points
if nargin > 5
    for k = 1:length(nanowires)
        % Particular nanowire
        nw = nanowires{k};
    
        % Load constants
        xc    = nw.xc;
        r     = nw.r;
        
        % Compute points to compute error in
        nw.phi = linspace(0,2*pi,n)';
        nw.x = r*cos(nw.phi) + xc;
        nw.y = r*sin(nw.phi) + r;
        nanowires{k} = nw;
    end
end

num_nanowires = length(nanowires);

% Make vectors ready for shared points
phi = [];
coord = struct;
coord.x = [];
coord.y = [];

% Get auxiliary sources and test points
for k = 1:num_nanowires
    % Particular nanowire
    nw = nanowires{k};
    
    % Combine all points for computation
    phi = [phi; nw.phi];
    coord.x = [coord.x; nw.x];
    coord.y = [coord.y; nw.y];
end

M = length(coord.x);

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario);
[Eref, Href] = reflected_fields(coord, scenario);

% Compute scattered fields
Escat = zeros(M,3);
Hscat = zeros(M,3);
if num_computations > 1
    coord_temp = struct;
    for k = 1:num_computations
        coord_temp.x = coord_int.x(k,:);
        coord_temp.y = coord_int.y(k,:);
        Escat(:,3) = Escat(:,3) + Ez_scat_matrix(coord, coord_temp) * C(:,k);
        Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(coord, coord_temp) * C(:,k);
        Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(coord, coord_temp) * C(:,k);
    end
else
    Escat(:,3) = Escat(:,3) + Ez_scat_matrix(coord, coord_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(coord, coord_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(coord, coord_int) * C;
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

if num_computations == 1

    Etot_inside = zeros(M,3);
    Htot_inside = zeros(M,3);

    Etot_inside(:,3) = Ez_tot_inside_matrix(coord, coord_ext) * D;
    Htot_inside(:,1) = Hx_tot_inside_matrix(coord, coord_ext) * D;
    Htot_inside(:,2) = Hy_tot_inside_matrix(coord, coord_ext) * D;

    E_error = abs(Etot(:,3) - Etot_inside(:,3))./abs(Etot_inside(:,3));

    H_left  = -sin(phi) .* Htot(:,1)        + cos(phi) .* Htot(:,2);
    H_right = -sin(phi) .* Htot_inside(:,1) + cos(phi) .* Htot_inside(:,2);

    H_error = abs(H_left - H_right)./abs(H_right);

    E_error = reshape(E_error,[],num_nanowires);
    H_error = reshape(H_error,[],num_nanowires);
else
    E_error = [];
    H_error = [];
    coord_temp = struct;
    for k = 1:num_nanowires
        % Particular nanowire
        nw = nanowires{k};
        x = nw.x;
        y = nw.y;
        phi = nw.phi;
        M = length(x);
        coord_temp.x = coord_ext.x(k,:);
        coord_temp.y = coord_ext.y(k,:);
        coord.x = x;
        coord.y = y;
    
        Etot_inside = zeros(M,3);
        Htot_inside = zeros(M,3);
    
        Etot_inside(:,3) = Ez_tot_inside_matrix(coord, coord_temp) * D(:,k);
        Htot_inside(:,1) = Hx_tot_inside_matrix(coord, coord_temp) * D(:,k);
        Htot_inside(:,2) = Hy_tot_inside_matrix(coord, coord_temp) * D(:,k);

        start = (k-1)*M + 1;
        stop  = k*M;

        error = abs(Etot(start:stop,3) - Etot_inside(:,3))./abs(Etot(start:stop,3));
    
        E_error = [E_error; error'];
    
        H_left  = -sin(phi) .* Htot(start:stop,1)        + cos(phi) .* Htot(start:stop,2);
        H_right = -sin(phi) .* Htot_inside(:,1) + cos(phi) .* Htot_inside(:,2);
        error = abs(H_left - H_right)./abs(H_left);

        H_error = [H_error; error'];
    end

    E_error = E_error';
    H_error = H_error';


end
