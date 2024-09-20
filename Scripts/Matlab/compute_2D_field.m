function [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_2D_field(X,Y,nanowires,C,D,x,y,phi,x_int,y_int,x_ext,y_ext)

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants_nanowires();

% Get size of field to compute
[n,m] = size(X);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
M = n*m;

% Get the number of seperate computations that have been done
[cl, num_computations] = size(C)

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

% Compute scattered fields
Escat = zeros(M,3);
Hscat = zeros(M,3);
if num_computations > 1
    for k = 1:num_computations
        Escat(:,3) = Escat(:,3) + Ez_scat_matrix(X,Y,x_int(k,:),y_int(k,:)) * C(:,k);
        Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(X,Y,x_int(k,:),y_int(k,:)) * C(:,k);
        Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(X,Y,x_int(k,:),y_int(k,:)) * C(:,k);
    end
else
    size(Ez_scat_matrix(X,Y,x_int,y_int))
    size(C)
    Escat(:,3) = Escat(:,3) + Ez_scat_matrix(X,Y,x_int,y_int) * C;
    Hscat(:,1) = Hscat(:,1) + Hx_scat_matrix(X,Y,x_int,y_int) * C;
    Hscat(:,2) = Hscat(:,2) + Hy_scat_matrix(X,Y,x_int,y_int) * C;
end

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;

% Find the interior of the nanowires
numerical = @(x,y) sqrt(x.^2+y.^2);
find_interior = zeros(m*n,1);
for j = 1:length(nanowires)
    nw = nanowires{j};
    dist = numerical(X - nw.xc, Y - nw.r);
    find_interior = find_interior + (dist < nw.r);
end

% Remove the fields form the interior of the nanowires
Etot(find_interior > 0,:) = 0;
Htot(find_interior > 0,:) = 0;


% Compute total fields inside nanowires
if num_computations > 1
    for k = 1:num_computations
        Etot_inside = zeros(M,3);
        Htot_inside = zeros(M,3);
        Etot_inside(:,3) = Ez_tot_inside_matrix(X,Y,x_ext(k,:),y_ext(k,:)) * D(:,k);
        Htot_inside(:,1) = Hx_tot_inside_matrix(X,Y,x_ext(k,:),y_ext(k,:)) * D(:,k);
        Htot_inside(:,2) = Hy_tot_inside_matrix(X,Y,x_ext(k,:),y_ext(k,:)) * D(:,k);
    
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

    Etot_inside(:,3) = Ez_tot_inside_matrix(X,Y,x_ext,y_ext) * D;
    Htot_inside(:,1) = Hx_tot_inside_matrix(X,Y,x_ext,y_ext) * D;
    Htot_inside(:,2) = Hy_tot_inside_matrix(X,Y,x_ext,y_ext) * D;

    % Remove field outside of nanowires
    Etot_inside(~find_interior,:) = 0;
    Htot_inside(~find_interior,:) = 0;
    
    % Add total field inside nanowires
    Etot = Etot + Etot_inside;
    Htot = Htot + Htot_inside;

end

Etot = reshape(Etot,n,m);
Htot = reshape(Htot,n,m);
Einc = reshape(Einc,n,m);
Hinc = reshape(Hinc,n,m);
Eref = reshape(Eref,n,m);
Href = reshape(Href,n,m);
Escat = reshape(Escat,n,m);
Hscat = reshape(Hscat,n,m);


return




% Do setup for easier computations
H02 = @(z) besselh(0,2,z);
H12 = @(z) besselh(1,2,z);


% Compute scattered fields from each nanowire
Escat = zeros(n,m,3);
Hscat = zeros(n,m,3);
wb = waitbar(0,'Computing scattered fields...');
br = 0;
for j = 1:length(nanowires)
    nw = nanowires{j};
    for k = 1:length(nw.C)
        abs_int      = numerical(X - nw.x_int(k), Y - nw.y_int(k));
        abs_int_ref  = numerical(X - nw.x_int(k), Y + nw.y_int(k));
        Escat(:,:,3) = Escat(:,:,3) + nw.C(k) * (H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref));
        val          =           1./abs_int     .* H12(k0 * abs_int);
        val_ref      = Gamma_r * 1./abs_int_ref .* H12(k0 * abs_int_ref);
        cons         = 1i/eta0 * nw.C(k);
        Hscat(:,:,1) = Hscat(:,:,1) - cons * (val .* (Y - nw.y_int(k)) + val_ref .* (Y + nw.y_int(k)));
        Hscat(:,:,2) = Hscat(:,:,2) + cons * (val .* (X - nw.x_int(k)) + val_ref .* (X - nw.x_int(k)));
        br=br+1;
        waitbar(br/(length(nanowires)*length(nw.C)),wb);
    end
end
close(wb);

% Get the total fields
Etot = Einc + Eref + Escat;
Htot = Hinc + Href + Hscat;


% Compute total field inside nanowires
wb = waitbar(0,'Computing total fields inside the nanowires...');
br = 0;
for j = 1:length(nanowires)
    nw = nanowires{j};
    Etot_inside = zeros(n,m,3);
    Htot_inside = zeros(n,m,3);
    if combined
        for s = 1:length(nanowires)
            nw_loc = nanowires{s};
            for k = 1:length(nw_loc.D)
                abs_ext            = numerical(X - nw_loc.x_ext(k), Y - nw_loc.y_ext(k));
                Etot_inside(:,:,3) = Etot_inside(:,:,3) + nw_loc.D(k) * H02(nw_loc.k1*abs_ext);
                cons = 1i * nw_loc.n1/eta0 * nw_loc.D(k) * 1./abs_ext;
                Hval = H12(nw_loc.k1*abs_ext);
                Htot_inside(:,:,1) = Htot_inside(:,:,1) - cons .* Hval .* (Y - nw_loc.y_ext(k));
                Htot_inside(:,:,2) = Htot_inside(:,:,2) + cons .* Hval .* (X - nw_loc.x_ext(k));
                br=br+1;
                waitbar(br/(length(nanowires)*length(nanowires)*length(nw_loc.D)),wb);
            end
        end

    else
        for k = 1:length(nw.D)
            abs_ext            = numerical(X - nw.x_ext(k), Y - nw.y_ext(k));
            Etot_inside(:,:,3) = Etot_inside(:,:,3) + nw.D(k) * H02(nw.k1*abs_ext);
            cons = 1i * nw.n1/eta0 * nw.D(k) * 1./abs_ext;
            Hval = H12(nw.k1*abs_ext);
            Htot_inside(:,:,1) = Htot_inside(:,:,1) - cons .* Hval .* (Y - nw.y_ext(k));
            Htot_inside(:,:,2) = Htot_inside(:,:,2) + cons .* Hval .* (X - nw.x_ext(k));
            br=br+1;
            waitbar(br/(length(nanowires)*length(nw.D)),wb);
        end
    end

    dist = numerical(X - nw.xc, Y - nw.r);

    % Remove field outside of nanowires
    for i = 1:3
        Etemp = Etot_inside(:,:,i);
        Htemp = Htot_inside(:,:,i);
        Etemp(dist >= nw.r) = 0;
        Htemp(dist >= nw.r) = 0;
        Etot_inside(:,:,i) = Etemp;
        Htot_inside(:,:,i) = Htemp;
    end

    % Add total field inside nanowires
    Etot = Etot + Etot_inside;
    Htot = Htot + Htot_inside;
end
close(wb);




