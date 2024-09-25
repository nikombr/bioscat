function [nanowires, X, Y] = setup_settings(choice, N, num_grid_points, view, dist)

Y = linspace(0,16,num_grid_points)*10^(-7);

if choice == 1
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire(N);
    X = linspace(-6,14,num_grid_points)*10^(-7);
elseif choice == 2
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire(N);
    nanowires{2} = setup_nanowire(N,7*10^(-7),1.2*nanowires{1}.r);
    X = linspace(-6,14,num_grid_points)*10^(-7);
elseif choice == 3
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire(N);
    nanowires{2} = setup_nanowire(N,2.45*nanowires{1}.r,1.2*nanowires{1}.r);
    X = linspace(-6,14,num_grid_points)*10^(-7);
elseif choice == 4
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire(N);
    nanowires{2} = setup_nanowire(N, 6*10^(-7),1.2*nanowires{1}.r);
    nanowires{3} = setup_nanowire(N,10*10^(-7),0.6*nanowires{1}.r);
    nanowires{4} = setup_nanowire(N,14*10^(-7),1.1*nanowires{1}.r);
    X = linspace(-6,20,num_grid_points)*10^(-7);
    Y = linspace(0,21,num_grid_points)*10^(-7);
elseif choice == 5
    nanowires = cell(1,1);
    nanowires{1} = setup_nanowire(N);
    nanowires{2} = setup_nanowire(N,2*nanowires{1}.r + dist,nanowires{1}.r);
    %X = linspace(-6,14,num_grid_points)*10^(-7);
    X = linspace(-6,46,num_grid_points)*10^(-7);
    Y = 0*X;
else
    disp("Please choose 1, 2, 3, 4 or 5.")
    return
end

if strcmp("far",view)
    Y = Y + 3*10^(-2);
end

for k = 1:length(nanowires)
    nanowires{k}.N = N;
end
