set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')
%% Compute results

clear; clc; close all; 

view = "close"; % close or far

settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
choice = 2;
num_grid_points = 100;
N = 100; % Number of auxiliary sources, either interior or exterior, or test points

Y = linspace(0,16,num_grid_points)*10^(-7);

if choice == 1
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire();
    X = linspace(-6,14,num_grid_points)*10^(-7);
elseif choice == 2
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire();
    nanowires{2} = setup_nanowire();
    nanowires{2}.xc = 7*10^(-7);
    nanowires{2}.r = 1.2*nanowires{2}.r;
    X = linspace(-6,14,num_grid_points)*10^(-7);

elseif choice == 3
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire();
    nanowires{2} = setup_nanowire();
    nanowires{2}.r = 1.2*nanowires{2}.r;
    nanowires{2}.xc = 1.1*(nanowires{1}.r + nanowires{2}.r);
    X = linspace(-6,14,num_grid_points)*10^(-7);

elseif choice == 4
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire();
    nanowires{2} = setup_nanowire();
    nanowires{2}.xc = 6*10^(-7);
    nanowires{2}.r = 1.2*nanowires{2}.r;
    nanowires{3} = setup_nanowire();
    nanowires{3}.xc = 10*10^(-7);
    nanowires{3}.r = 0.7*nanowires{2}.r;
    nanowires{4} = setup_nanowire();
    nanowires{4}.xc = 14*10^(-7);
    nanowires{4}.r = 1.1*nanowires{2}.r;
    X = linspace(-6,20,num_grid_points)*10^(-7);
    Y = linspace(0,21,num_grid_points)*10^(-7);

else
    disp("Please choose 1, 2, 3 or 4.")
end

if strcmp("far",view)
    Y = Y + 3*10^(-2);
end


for k = 1:length(nanowires)
    nanowires{k}.N = N;
end

C = [];
D = [];
x_int = [];
y_int = [];
x_ext = [];
y_ext = [];
for k = 1:length(nanowires)
    temp = {nanowires{k}};
    [c,d,x,y,phi,t_x_int,t_y_int,t_x_ext,t_y_ext] = forward_nanowires(temp,N);
    C = [C c];
    D = [D d];
    x_int = [x_int; t_x_int];
    y_int = [y_int; t_y_int];
    x_ext = [x_ext; t_x_ext];
    y_ext = [y_ext; t_y_ext];
end
size(x_int)
size(C)
[Xmesh,Ymesh] = meshgrid(X,Y);
[Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_2D_field(Xmesh,Ymesh,nanowires,C,D,x,y,phi,x_int,y_int,x_ext,y_ext);

filename = sprintf('../../Results/nanowires/seperate_computation_%s_%s_N_%d.mat',view,settings{choice},N);
save(filename,'Etot','Htot','Einc','Hinc','Eref','Href','Escat','Hscat','X','Y','nanowires');

[C,D,x,y,phi,x_int,y_int,x_ext,y_ext] = forward_nanowires(nanowires,N);
[Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_2D_field(Xmesh,Ymesh,nanowires,C,D,x,y,phi,x_int,y_int,x_ext,y_ext);

filename = sprintf('../../Results/nanowires/combined_computation_%s_%s_N_%d.mat',view,settings{choice},N);
save(filename,'Etot','Htot','Einc','Hinc','Eref','Href','Escat','Hscat','X','Y','nanowires');


%%

clear; close all; clc;

settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
choice = 2;


filename = sprintf('../../Results/nanowires/seperate_computation_close_%s_N_100.mat',settings{choice});
load(filename);
combined = false;
[E_error, H_error] = compute_transmission_error(nanowires,combined);


plot(E_error')

