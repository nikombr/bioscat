function [b1, b2, a1, a2, a3, a4] = forward_section(coord, coord_int, coord_ext, scenario, lambda0, section, n_x, n_y)

% When solving the linear system we are in the near field, so we cannot do far field approximations.
far_field_approximation = false;

% Compute incident and reflected fields
[Einc, Hinc] =  incident_fields(coord, scenario, lambda0);
[Eref, Href] = reflected_fields(coord, scenario, lambda0);

% Get incident fields
Ez_inc = Einc(:,3);
Hx_inc = Hinc(:,1);

% Get reflected fields
%Ez_ref = Eref(:,3);
%Hx_ref = Href(:,1);

% Compute scattered fields
%Ez_scat = Ez_scat_matrix(coord, coord_int);
%Hx_scat = Hx_scat_matrix(coord, coord_int);
%Hy_scat = Hy_scat_matrix(coord, coord_int);

% Compute total fields inside the nanostructure
%Ez_tot_inside = Ez_tot_inside_matrix(coord, coord_ext);
%Hx_tot_inside = Hx_tot_inside_matrix(coord, coord_ext);
%Hy_tot_inside = Hy_tot_inside_matrix(coord, coord_ext);

[Escat, Hscat] = scattered_fields(coord, coord_int, scenario, lambda0, far_field_approximation);

if strcmp(section, "top")

    % Setup vector
    b1 = - Ez_inc - Ez_ref;
    b2 = n_y .* (Hx_inc + Hx_ref);  % z component

    % Setup matrix
    a1 = Ez_scat;
    a2 = -Ez_tot_inside;
    a3 =   n_x .* Hy_scat - n_y .* Hx_scat;
    a4 = - n_x .* Hy_tot_inside + n_y .* Hx_tot_inside;

elseif strcmp(section, "right")

    % Setup vector
    b1 = - Ez_inc - Ez_ref;
    b2 = 0*Ez_inc;

    % Setup matrix
    a1 = Ez_scat;
    a2 = -Ez_tot_inside;
    a3 = Hy_scat;
    a4 = -Hy_tot_inside;

elseif strcmp(section, "bottom")

    % Setup vector
    b1 = - Ez_inc - Ez_ref;
    b2 = Hx_inc + Hx_ref;

    % Setup matrix
    a1 = Ez_scat;
    a2 = -Ez_tot_inside;
    a3 = Hx_scat;
    a4 = -Hx_tot_inside;

elseif strcmp(section, "left")
    
    % Setup vector
    b1 = - Ez_inc - Ez_ref;
    b2 = 0*Ez_inc;

    % Setup matrix
    a1 = Ez_scat;
    a2 = -Ez_tot_inside;
    a3 = Hy_scat;
    a4 = -Hy_tot_inside;

else
    fprintf("Please inpute 'top', 'right', 'bottom' or 'left'. You wrote '%s'.\n",section);
end