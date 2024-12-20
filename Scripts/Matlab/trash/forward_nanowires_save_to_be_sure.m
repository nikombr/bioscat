function [nanowires] = forward_nanowires(nanowires)

m = length(nanowires);

% Load general constants
[eta0, n0, ns, lambda0, Gamma_r, Gamma_t, k0, ks] = load_constants_nanowires();

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
    N     = nw.N;
    xc    = nw.xc;
    r     = nw.r;
    k1    = nw.k1;
    alpha = nw.alpha;
    n1    = nw.n1;
    
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
    %phi = [phi; nw.phi];
    %x = [x; nw.x];
    %y = [y; nw.y];
    x_int = [x_int nw.x_int];
    y_int = [y_int nw.y_int];
    x_ext = [x_ext nw.x_ext];
    y_ext = [y_ext nw.y_ext];
    
end

% Make vector and matrix vector for shared system
b = [];
A = [];

for k = 1:m
    % Particular nanowire
    nw = nanowires{k};

    % Load constants
    N     = nw.N;
    k1    = nw.k1;
    n1    = nw.n1;
    phi   = nw.phi;
    x     = nw.x;
    y     = nw.y;
    
    % Setup vector
    b1 = - exp(1i*k0*y) - Gamma_r * exp(-1i*k0*y);
    b2 = sin(phi)/eta0 .* (-exp(1i*k0*y) + Gamma_r * exp(-1i*k0*y));
    b = [b; b1; b2];
    
    % Do precomputations
    H02 = @(z) besselh(0,2,z);
    H12 = @(z) besselh(1,2,z);
    
    numerical = @(x,y) sqrt(x.^2+y.^2);
    
    abs_int     = numerical(x - x_int, y - y_int);
    abs_int_ref = numerical(x - x_int, y + y_int);
    abs_ext     = numerical(x - x_ext, y - y_ext);
    
    % Steup matrix
    a1 = H02(k0*abs_int) + Gamma_r*H02(k0*abs_int_ref);
    a2 = -H02(k1*abs_ext);
    a3 = 1i/eta0 * (1./abs_int     .* H12(k0*abs_int)      .* (sin(phi).*(y - y_int)+cos(phi).*(nw.x - x_int)) + ...
                    1./abs_int_ref .* H12(k0*abs_int_ref)  .* (sin(phi).*(y + y_int)+cos(phi).*(x - x_int)));
    a4 = -1i*n1/eta0 * 1./abs_ext  .* (sin(phi).*(y - y_ext) + cos(phi).*(x - x_ext)).*H12(k1*abs_ext);
    temp = [a1 a2; a3 a4];
    A = [A; temp];
    

end

% Solve linear system
c = A\b;
length(c)
% Look at the error of the linear system
error = max(abs(A*c - b))

% Save result for each nanowire
for k = 1:m
    Cstart  = (k-1)*N+1
    Cend    = (k-1)*N+N
    Dstart  = (k-1)*N + m*N + 1
    Dend    = (k-1)*N + m*N + N
    nanowires{k}.C = c(Cstart:Cend);
    nanowires{k}.D = c(Dstart:Dend);
end


