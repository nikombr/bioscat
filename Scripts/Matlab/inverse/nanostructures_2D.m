
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

setup = "2D";
total_grid_points = 1000;
protein_structure = 'demoleus2x2';  % Retinin2x2 demoleus2x2
delta = 0.01;
num_segments = 10;
data_quality = "clean";

covfunc = "matern";
arg1 = "4";
arg2 = "0.8";
num = 7000;

preConditionedCrankNicholsonScheme(setup,num,protein_structure,total_grid_points,delta,covfunc,arg1,arg2,num_segments,data_quality)

%delete(gcp('nocreate'))