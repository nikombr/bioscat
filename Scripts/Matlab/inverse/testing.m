
%parpool;

%%

clear; close all; clc;

dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
%dir = '/zhome/00/b/147112/bioscat/';

addpath(sprintf('%s/Scripts/Matlab/utils/',dir))
addpath(sprintf('%s/Scripts/Matlab/forward/nanostructures_2D',dir))
addpath(sprintf('%s/Scripts/Matlab/forward/fields_2D',dir))

warning('off', 'MATLAB:singularMatrix');  % Suppresses singular matrix warnings
warning('off', 'MATLAB:nearlySingularMatrix')

setup = "test";
total_grid_points = 1000;
protein_structure = 'demoleus2x2';  % Retinin2x2 demoleus2x2
delta = 0.01;

covfunc = "matern";
deltas = [0.001 0.01 0.1 0.2 0.3];
arg1 = {"1","2","4"};
arg2 = {"0.1","0.2","0.4","0.8","2"};

for delta = deltas
    for i = 1:length(arg1)
        for j = 1:length(arg2)
            close all;
            clc;
            preConditionedCrankNicholsonScheme(setup,protein_structure,total_grid_points,delta,covfunc,arg1{i},arg2{j})
        end
    end
end

%delete(gcp('nocreate'))