function [Eref, Href] = reflected_fields(coord, scenario, lambda0)

y = coord.y;

M = length(y);

% Allocation
Eref = zeros(M,3);
Href = zeros(M,3);

% Load constants
const   = load_constants(lambda0);
k0      = const.k0;
eta0    = const.eta0;
Gamma_r = const.Gamma_r;

% Compute fields
if scenario == 1
    
    Eref(:,3) =          Gamma_r * exp(-1i*k0*y);
    Href(:,1) = 1/eta0 * Gamma_r * exp(-1i*k0*y);

    1/eta0*Gamma_r * exp(0)
    1/eta0*Gamma_r * cos(0)

elseif scenario == 2

    Eref(:,1) =            Gamma_r * exp(-1i*k0*y);
    Href(:,3) = - 1/eta0 * Gamma_r * exp(-1i*k0*y);

else
    fprintf("You have to input 1 or 2 for scenario, but you wrote %d\n",scenario)
end