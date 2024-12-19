clear; close all; clc; 


% Load modules
import com.comsol.model.*
import com.comsol.model.util.*



% Set variables
%lambda=325; % nm
%beta = 0;
%protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2
protein_structures = {'Retinin2x2','demoleus2x2'};
for p = 1:2
    protein_structure = protein_structures{p}
    for beta = [0, 30, 60, 90]
        beta
        for lambda = [325 425 525 625]
            lambda
            generateComsolData(protein_structure, beta, lambda)
        end
    end
end

function generateComsolData(protein_structure, beta, lambda, obs_grid)

% Set default values
if nargin < 4
    obs_grid = 500;
end

% Directory data and results
dir = "../../../../../../../work3/s194146/bioscatdata";

% Load Comsol model
model = mphload(sprintf('%s/../forward61_%s.mph',dir,protein_structure));

% Set the parameter in the COMSOL model
param_name = 'lambda_0';    % Name of the parameter in the COMSOL model
param_value = lambda*1e-9;          % New value for the parameter
model.param.set(param_name, param_value);
param_name = 'beta';    % Name of the parameter in the COMSOL model
param_value = beta*pi/180;          % New value for the parameter
model.param.set(param_name, param_value);

% Solve
model.sol('sol1').runAll;

% Define angles to find the far field in
phi = linspace(0,pi,obs_grid);

% Setup points to evaluate comsol model in 
dc = 0.001;
y = [0:dc:4 (-2.5:dc:2.5)*0+4 4:-dc:0]*1e-6;
x = [(0:dc:4)*0-2.5 -2.5:dc:2.5 (4:-dc:0)*0+2.5]*1e-6;
coord_comsol = [x; y];

% Get angles we find the solution in from comsol
phi_comsol = atan2(y,x);

% Get result from comsol
result = mphinterp(model, 'emw.normEfar', 'coord', coord_comsol);

% Rempve repeated points
[phi_comsol,i,j] = unique(phi_comsol);
result = result(i);

% Get solution for specified angles
val = interp1(phi_comsol, result, phi);

filename = sprintf("%s/Data/comsol/%s/far_field_pattern/radiation_pattern_beta_%d_lambda_%d_obs_grid_%d.mat",dir,protein_structure,beta,lambda,obs_grid);
save(filename,'phi','val')

end