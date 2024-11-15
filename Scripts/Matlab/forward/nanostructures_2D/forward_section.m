function [b1, b2, a1, a2, a3, a4] = forward_section(coord, coord_int, coord_ext, scenario, lambda0, section, n_x, n_y)

% When solving the linear system we are in the near field, so we cannot do far field approximations.
far_field_approximation = false;

% Compute fields
[Einc,        Hinc]        =  incident_fields(coord,            scenario, lambda0);
[Eref,        Href]        = reflected_fields(coord,            scenario, lambda0);
[Escat,       Hscat]       = scattered_fields(coord, coord_int, scenario, lambda0, far_field_approximation);
[Etot_inside, Htot_inside] =  interior_fields(coord, coord_ext, scenario, lambda0);
Ez = Etot_inside(:,1,3)
Hx = Htot_inside(:,1,1)
Hy = Htot_inside(:,1,2)

if scenario == 2
    % Swap fields in scenario 2 as the linear system can be set up in the
    % exact same way as scenario with the exception that the fields are
    % reveresed.
    [Hinc,        Einc]        = deal(Einc,        Hinc);
    [Href,        Eref]        = deal(Eref,        Href);
    [Hscat,       Escat]       = deal(Escat,       Hscat);
    [Htot_inside, Etot_inside] = deal(Etot_inside, Htot_inside);
end

% Get incident fields
%Ez_inc = Einc(:,3);
%Hx_inc = Hinc(:,1);

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

b1 = - Einc(:,3) - Eref(:,3);
a1 = Escat(:,:,3);
a2 = -Etot_inside(:,:,3);

if strcmp(section, "top")

    b2 =   n_y .* (Hinc(:,1) + Href(:,1));  % z component
    a3 =   n_x .* Hscat(:,:,2) - n_y .* Hscat(:,:,1);
    a4 = - n_x .* Htot_inside(:,:,2) + n_y .* Htot_inside(:,:,1);

elseif strcmp(section, "right") || strcmp(section, "left") 

    b2 = 0*Hinc(:,1);
    a3 = Hscat(:,:,2);
    a4 = -Htot_inside(:,:,2);

elseif strcmp(section, "bottom")

    b2 = Hinc(:,1) + Href(:,1);
    a3 = -Hscat(:,:,1);
    a4 = Htot_inside(:,:,1);

else
    fprintf("Please inpute 'top', 'right', 'bottom' or 'left'. You wrote '%s'.\n",section);
end



% 
% if strcmp(section, "top")
% 
%     % Setup vector
%     b1 = - Ez_inc - Ez_ref;
%     b2 = n_y .* (Hx_inc + Hx_ref);  % z component
% 
%     % Setup matrix
%     a1 = Ez_scat;
%     a2 = -Ez_tot_inside;
%     a3 =   n_x .* Hy_scat - n_y .* Hx_scat;
%     a4 = - n_x .* Hy_tot_inside + n_y .* Hx_tot_inside;
% 
% elseif strcmp(section, "right")
% 
%     % Setup vector
%     b1 = - Ez_inc - Ez_ref;
%     b2 = 0*Ez_inc;
% 
%     % Setup matrix
%     a1 = Ez_scat;
%     a2 = -Ez_tot_inside;
%     a3 = Hy_scat;
%     a4 = -Hy_tot_inside;
% 
% elseif strcmp(section, "bottom")
% 
%     % Setup vector
%     b1 = - Ez_inc - Ez_ref;
%     b2 = Hx_inc + Hx_ref;
% 
%     % Setup matrix
%     a1 = Ez_scat;
%     a2 = -Ez_tot_inside;
%     a3 = Hx_scat;
%     a4 = -Hx_tot_inside;
% 
% elseif strcmp(section, "left")
% 
%     % Setup vector
%     b1 = - Ez_inc - Ez_ref;
%     b2 = 0*Ez_inc;
% 
%     % Setup matrix
%     a1 = Ez_scat;
%     a2 = -Ez_tot_inside;
%     a3 = Hy_scat;
%     a4 = -Hy_tot_inside;
% 
% else
%     fprintf("Please inpute 'top', 'right', 'bottom' or 'left'. You wrote '%s'.\n",section);
% end