
%parpool;

%%

clear; close all; clc;

dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
dir = '/zhome/00/b/147112/bioscat/';

addpath(sprintf('%s/Scripts/Matlab/utils/',dir))
addpath(sprintf('%s/Scripts/Matlab/forward/nanostructures_2D',dir))
addpath(sprintf('%s/Scripts/Matlab/forward/fields_2D',dir))

warning('off', 'MATLAB:singularMatrix');  % Suppresses singular matrix warnings
warning('off', 'MATLAB:nearlySingularMatrix')

preConditionedCrankNicholsonScheme(@testSetup)

%delete(gcp('nocreate'))