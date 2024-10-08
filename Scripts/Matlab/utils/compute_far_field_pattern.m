function [phi, PE, PH] = compute_far_field_pattern(dim,num,r,varargin)

if r < 10^(-3)
    print("I don't think this qualifies as the far field :(");
    return
end

% Define norm
normfunc = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

if dim == 2 % 2D computations
    
    phi = linspace(0,pi,num);
    
    X = r*cos(phi);
    Y = r*sin(phi);

    coord = struct;
    coord.x = X;
    coord.y = Y;

    [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(coord,varargin{:});


    % Compute far-field pattern
    PE = normfunc(abs(Escat));
    PH = normfunc(abs(Hscat));
    phi = phi';


elseif dim == 3







end