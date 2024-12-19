%clear; close all; clc; 


% Load modules
import com.comsol.model.*
import com.comsol.model.util.*

% Set variables
%lambda=325; % nm
%beta = 0;
%protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2
protein_structures = {'Retinin2x2','demoleus2x2'};
idx = 2;

%fine_part = 5
%protein_structure = protein_structures{idx}
%phi = pi/2;
%betas = 0:10:90;
%if fine_part == 1
%    lambdas = 250:5:750;
%elseif fine_part == 2
%    lambdas = 251:5:750;
%elseif fine_part == 3
%    lambdas = 252:5:750;
%elseif fine_part == 4
%    lambdas = 253:5:750;
%elseif fine_part == 5
%    lambdas = 254:5:750;
%end
%folder = sprintf("one_observation_ultra_fine_%d",fine_part);

%generateComsolData(protein_structure, phi, betas, lambdas, folder);

%protein_structure = protein_structures{idx}
%phi = pi/2;
%betas = 0:10:90;
%lambdas = 250:1:750;
%folder = "one_observation_ultra_fine";

%generateComsolData(protein_structure, phi, betas, lambdas, folder);

%protein_structure = protein_structures{idx};
%phi = pi/2;
%betas = 0:10:90;
%lambdas = 250:5:750;
%folder = "one_observation_fine";

%generateComsolData(protein_structure, phi, betas, lambdas, folder);

protein_structure = protein_structures{idx}
phi = linspace(0,pi,200);
betas = 0:10:90;
lambdas = 250:750;
lambdas = [325];
folder = "angle_resolved_one_wavelength";

generateComsolData(protein_structure, phi, betas, lambdas, folder);

%protein_structure = protein_structures{idx}
%phi = linspace(0,pi,200);
%betas = 0:10:90;
%lambdas = 250:750;
%lambdas = [325 425 525 625];
%folder = "angle_resolved";

%generateComsolData(protein_structure, phi, betas, lambdas, folder);

%protein_structure = protein_structures{idx}
%phi = pi/2;
%betas = 0:10:90;
%lambdas = 250:10:750;
%folder = "one_observation";

%generateComsolData(protein_structure, phi, betas, lambdas, folder);

quit

function generateComsolData(protein_structure, phi, betas, lambdas, folder)

% Directory data and results
dir = "/work3/s194146/bioscatdata";

% Load Comsol model
model = mphload(sprintf('%s/../forward61_%s.mph',dir,protein_structure));
reflectance = zeros(length(lambdas), length(betas), length(phi));
count = 0;
for l = 1:length(lambdas)
    lambda = lambdas(l);
    for b = 1:length(betas)
        beta = betas(b);
        fprintf("\t%d/%d\r",count,length(lambdas)*length(betas))
        tic;
        % Set the parameter in the COMSOL model
        param_name = 'lambda_0';    % Name of the parameter in the COMSOL model
        param_value = lambda*1e-9;          % New value for the parameter
        model.param.set(param_name, param_value);
        param_name = 'beta';    % Name of the parameter in the COMSOL model
        param_value = beta*pi/180;          % New value for the parameter
        model.param.set(param_name, param_value);

        % Solve
        model.sol('sol1').runAll;

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
        reflectance(l,b,:) = interp1(phi_comsol, result, phi);
        toc;
        count = count + 1;
    end
end

filename = sprintf("%s/Data/artificial_data/comsol/%s/%s/reflectance.txt",dir,folder,protein_structure);
writematrix(reshape(reflectance,[],1), filename)
filename = sprintf("%s/Data/artificial_data/comsol/%s/%s/betas.txt",dir,folder,protein_structure);
writematrix(betas'*pi/180, filename)
filename = sprintf("%s/Data/artificial_data/comsol/%s/%s/lambdas.txt",dir,folder,protein_structure);
writematrix(lambdas'*1e-9, filename)
filename = sprintf("%s/Data/artificial_data/comsol/%s/%s/phi.txt",dir,folder,protein_structure);
writematrix(phi', filename)

end