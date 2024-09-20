function [C,D,x,y,phi,x_int,y_int,x_ext,y_ext] = forward_nanowires(nanowires,N)

addpath('compute_field_vectors')

m = length(nanowires);

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks, n1, alpha, k1] = load_constants_nanowires();

% Make vectors ready for shared points
phi = [];
x = [];
y = [];
x_int = [];
y_int = [];
x_ext = [];
y_ext = [];

% Get auxiliary sources and test points
for k = 1:m
    % Particular nanowire
    nw = nanowires{k};

    % Load constants
    xc    = nw.xc;
    r     = nw.r;
    
    % Compute auxiliary sources and test points
    nw.phi = linspace(0,2*pi,N)';
    nw.x = r*cos(nw.phi) + xc;
    nw.y = r*sin(nw.phi) + r;
    nw.x_int = alpha*r*cos(nw.phi') + xc;
    nw.y_int = alpha*r*sin(nw.phi') + r;
    nw.x_ext = 1/alpha*r*cos(nw.phi') + xc;
    nw.y_ext = 1/alpha*r*sin(nw.phi') + r;
    nanowires{k} = nw;

    % Combine all points for computation
    phi = [phi; nw.phi];
    x = [x; nw.x];
    y = [y; nw.y];
    x_int = [x_int nw.x_int];
    y_int = [y_int nw.y_int];
    x_ext = [x_ext nw.x_ext];
    y_ext = [y_ext nw.y_ext];
end

% Compute incident fields
Ez_inc = Ez_inc_vector(x,y);
Hx_inc = Hx_inc_vector(x,y);

% Compute reflected fields
Ez_ref = Ez_ref_vector(x,y);
Hx_ref = Hx_ref_vector(x,y);

    
% Setup vector
b1 = - Ez_inc - Ez_ref;
b2 = - sin(phi) .* (Hx_inc + Hx_ref);
%b1 = - exp(1i*k0*y) - Gamma_r * exp(-1i*k0*y);
%b2 = sin(phi)/eta0 .* (-exp(1i*k0*y) + Gamma_r * exp(-1i*k0*y));
b = [b1; b2];
    
% Do precomputations
%H02 = @(z) besselh(0,2,z);
%H12 = @(z) besselh(1,2,z);
    
%numerical = @(x,y) sqrt(x.^2+y.^2);
    
%abs_int     = numerical(x - x_int, y - y_int);
%abs_int_ref = numerical(x - x_int, y + y_int);
%abs_ext     = numerical(x - x_ext, y - y_ext);

% Compute scattered fields
Ez_scat = Ez_scat_matrix(x,y,x_int,y_int);
Hx_scat = Hx_scat_matrix(x,y,x_int,y_int);
Hy_scat = Hy_scat_matrix(x,y,x_int,y_int);

% Compute total fields inside the nanowires
Ez_tot_inside = Ez_tot_inside_matrix(x,y,x_ext,y_ext);
Hx_tot_inside = Hx_tot_inside_matrix(x,y,x_ext,y_ext);
Hy_tot_inside = Hy_tot_inside_matrix(x,y,x_ext,y_ext);
    
% Steup matrix
a1 = Ez_scat;
a2 = Ez_tot_inside;
a3 = -sin(phi) .* Hx_scat + cos(phi) .* Hy_scat;
a4 = -sin(phi) .* Hx_tot_inside + cos(phi) .* Hy_tot_inside;
%a1 = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);
%a2 = -H02(k1*abs_ext);
%a3 = 1i/eta0 * (1./abs_int     .* H12(k0*abs_int)      .* (sin(phi).*(y - y_int)+cos(phi).*(nw.x - x_int)) + ...
                %1./abs_int_ref .* H12(k0*abs_int_ref)  .* (sin(phi).*(y + y_int)+cos(phi).*(x - x_int)));
%a4 = -1i*n1/eta0 * 1./abs_ext  .* (sin(phi).*(y - y_ext) + cos(phi).*(x - x_ext)).*H12(k1*abs_ext);
A = [a1 a2; a3 a4];

% Solve linear system
c = A\b;

% Save result
C = c(1:length(c)/2);
D = c(length(c)/2+1:end);

% Look at the error of the linear system
error = max(abs(A*c - b));
fprintf("\nThe error from solving the linear system was %.4e\n\n",error);


