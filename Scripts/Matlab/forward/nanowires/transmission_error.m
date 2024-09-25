function [E_error, H_error] = transmission_error(C,D,x_int,y_int,x_ext,y_ext,nanowires,n)

% Get the number of seperate computations that have been done
[~, num_computations] = size(C);

% n is number of points uniformly distributed over the circles, if it is
% not provided, the error is computed in the test points
if nargin > 7
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

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants();


% nake vectors ready for shared points
phi = [];
x = [];
y = [];

% Get auxiliary sources and test points
for k = 1:num_nanowires
    % Particular nanowire
    nw = nanowires{k};
    
    % Combine all points for computation
    phi = [phi; nw.phi];
    x = [x; nw.x];
    y = [y; nw.y];
end

M = length(x);

% Compute incident fields
Einc = zeros(M,3);
Einc(:,3) = Ez_inc_vector(x,y);
Hinc = zeros(M,3);
Hinc(:,1) = Hx_inc_vector(x,y);

% Compute reflected field
Eref = zeros(M,3);
Eref(:,3) = Ez_ref_vector(x,y);
Href = zeros(M,3);
Href(:,1) = Hx_ref_vector(x,y);

% Compute scattered fields
Escat = zeros(M,3);
Hscat = zeros(M,3);
if num_computations > 1
    for k = 1:num_nanowires
        Escat(:,3) = Escat(:,3) + Ez_scat_matrix(x,y,x_int(k,:),y_int(k,:)) * C(:,k);
        Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(x,y,x_int(k,:),y_int(k,:)) * C(:,k);
        Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(x,y,x_int(k,:),y_int(k,:)) * C(:,k);
    end
else
    size(Ez_scat_matrix(x,y,x_int,y_int))
    size(C)
    Escat(:,3) = Escat(:,3) + Ez_scat_matrix(x,y,x_int,y_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(x,y,x_int,y_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(x,y,x_int,y_int) * C;
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

if num_computations == 1

    Etot_inside = zeros(M,3);
    Htot_inside = zeros(M,3);

    Etot_inside(:,3) = Ez_tot_inside_matrix(x,y,x_ext,y_ext) * D;
    Htot_inside(:,1) = Hx_tot_inside_matrix(x,y,x_ext,y_ext) * D;
    Htot_inside(:,2) = Hy_tot_inside_matrix(x,y,x_ext,y_ext) * D;

    E_error = abs(Etot(:,3) - Etot_inside(:,3));

    H_left  = -sin(phi) .* Htot(:,1)        + cos(phi) .* Htot(:,2);
    H_right = -sin(phi) .* Htot_inside(:,1) + cos(phi) .* Htot_inside(:,2);

    H_error = abs(H_left - H_right);

    E_error = reshape(E_error,[],num_nanowires);
    H_error = reshape(H_error,[],num_nanowires);
else
    E_error = [];
    H_error = [];
    for k = 1:num_nanowires
        % Particular nanowire
        nw = nanowires{k};
        x = nw.x;
        y = nw.y;
        phi = nw.phi;
        M = length(x);
    
        Etot_inside = zeros(M,3);
        Htot_inside = zeros(M,3);
    
        Etot_inside(:,3) = Ez_tot_inside_matrix(x,y,x_ext(k,:),y_ext(k,:)) * D(:,k);
        Htot_inside(:,1) = Hx_tot_inside_matrix(x,y,x_ext(k,:),y_ext(k,:)) * D(:,k);
        Htot_inside(:,2) = Hy_tot_inside_matrix(x,y,x_ext(k,:),y_ext(k,:)) * D(:,k);

        start = (k-1)*M + 1;
        stop  = k*M;

        error = abs(Etot(start:stop,3) - Etot_inside(:,3));
        size(error)
        size(E_error)
    
        E_error = [E_error; error'];
    
        H_left  = -sin(phi) .* Htot(start:stop,1)        + cos(phi) .* Htot(start:stop,2);
        H_right = -sin(phi) .* Htot_inside(:,1) + cos(phi) .* Htot_inside(:,2);
        error = abs(H_left - H_right);

        H_error = [H_error; error'];
    end

    E_error = E_error';
    H_error = H_error';


end


