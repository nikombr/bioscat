function nw = forward_one_nanowire(nw)

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks] = load_constants_nanowires();

% Load constants
N     = nw.N;
xc    = nw.xc;
r     = nw.r;
k1    = nw.k1;
alpha = nw.alpha;
n1    = nw.n1;

% Compute auxiliary sources and test points
phi = linspace(0,2*pi,N)';
x = r*cos(phi) + xc;
y = r*sin(phi) + r;
x_int = alpha*r*cos(phi') + xc;
y_int = alpha*r*sin(phi') + r;
x_ext = 1/alpha*r*cos(phi') + xc;
y_ext = 1/alpha*r*sin(phi') + r;

% Save for later
nw.phi   = phi;
nw.x     = x;
nw.y     = y;
nw.x_int = x_int;
nw.y_int = y_int;
nw.x_ext = x_ext;
nw.y_ext = y_ext;

% Setup vector
b1 = - exp(1i*k0*y) - Gamma_r * exp(-1i*k0*y);
b2 = sin(phi)/eta0 .* (-exp(1i*k0*y) + Gamma_r * exp(-1i*k0*y));
b = [b1; b2];

% Do precomputations
H02 = @(z) besselh(0,2,z);
H12 = @(z) besselh(1,2,z);

numerical = @(x,y) sqrt(x.^2+y.^2);

abs_int     = numerical(x-x_int,y-y_int);
abs_int_ref = numerical(x-x_int,y+y_int);
abs_ext     = numerical(x-x_ext,y-y_ext);

% Steup matrix
a1 = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);
a2 = -H02(k1*abs_ext);
a3 = 1i/eta0 * (1./abs_int     .* H12(k0*abs_int)     .* (sin(phi).*(y - y_int)+cos(phi).*(x - x_int)) + ...
                1./abs_int_ref .* H12(k0*abs_int_ref) .* (sin(phi).*(y + y_int)+cos(phi).*(x - x_int)));
a4 = -1i*n1/eta0 * 1./abs_ext .* (sin(phi).*(y - y_ext) + cos(phi).*(x - x_ext)).*H12(k1*abs_ext);
A = [a1 a2; a3 a4];

% Solve linear system
c = A\b;

% Look at the error of the linear system
error = max(abs(A*c - b))

% Save result
nw.C = c(1:N);
nw.D = c(N+1:end);

