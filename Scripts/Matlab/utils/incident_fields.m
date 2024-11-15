function [Einc, Hinc] = incident_fields(coord,scenario,varargin)

y = coord.y;

M = length(y);

% Allocation
Einc = zeros(M,3);
Hinc = zeros(M,3);

% Load constants
const = load_constants(varargin{:});
k0    = const.k0;
eta0  = const.eta0;

% Compute fields
if scenario == 1
    
    Einc(:,3) =            exp(1i*k0*y);
    Hinc(:,1) = - 1/eta0 * exp(1i*k0*y);

    - 1/eta0 * exp(0)
    - 1/eta0 * cos(0)


elseif scenario == 2

    Einc(:,1) =          exp(1i*k0*y);
    Hinc(:,3) = 1/eta0 * exp(1i*k0*y);


else
    fprintf("You have to input 1 or 2 for scenario, but you wrote %d\n",scenario)
end