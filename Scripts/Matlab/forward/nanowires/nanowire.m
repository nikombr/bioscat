set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')
%% Compute results in the near field

clear; clc; close all; 

addpath('../fields_2D')

num_grid_points = 300; % Number of grid points for plotting
N               = 100; % Number of auxiliary sources, either interior or exterior, or test points
views           = {"close", "far"};
settings        = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
far_field_approximation = false;

wb = waitbar(0,'Computing fields...');
br = 0;
for view_choice = 1:3
    if view_choice == 3
        far_field_approximation = true;
        view = views{2};
    else
        view = views{view_choice};
    end
    
    for choice = 1:4
        
        % Load settings for nanowires
        [nanowires, X, Y] = setup_settings(choice, N, num_grid_points, view);

        for computation = ["seperate", "combined"]
            
            if strcmp(computation,"seperate")

                % Run seperate computation
                tic;
                C = []; D = []; x_int = []; y_int = []; x_ext = []; y_ext = [];
                for k = 1:length(nanowires)
                    temp = {nanowires{k}};
                    [c,d,t_x_int,t_y_int,t_x_ext,t_y_ext] = forward(temp);
                    C = [C c];
                    D = [D d];
                    x_int = [x_int; t_x_int];
                    y_int = [y_int; t_y_int];
                    x_ext = [x_ext; t_x_ext];
                    y_ext = [y_ext; t_y_ext];
                end
                stop = toc;
                fprintf("\nIt took %.4f seconds to solve all the seperate linear systems.\n\n",stop)

            elseif strcmp(computation,"combined")

                % Run combined computation
                tic;
                [C,D,x_int,y_int,x_ext,y_ext] = forward(nanowires);
                stop = toc;
                fprintf("\nIt took %.4f seconds to solve the combined linear system.\n\n",stop)

            end
            
            % Save results
            extra = "";
            if strcmp(view,"far") && far_field_approximation
                extra = "_approximation";
            end
            filename = sprintf('../../../../Results/nanowires/%s_computation_%s%s_%s_N_%d.mat',computation,view,extra,settings{choice},N);
            [Xmesh,Ymesh] = meshgrid(X,Y);
            if strcmp(view,"far") && far_field_approximation
                tic;
                [Etot, Htot, Escat, Hscat, Einc, Hinc, Eref, Href] = compute_far_fields(Xmesh,Ymesh,C,x_int,y_int);
                stop = toc;
                fprintf("\nIt took %.4f seconds to compute the far fields.\n\n",stop)
                save(filename,'Etot','Htot','Einc','Hinc','Eref','Href','Escat','Hscat','X','Y','nanowires','C','D','x_int','y_int','x_ext','y_ext');
            else
                tic;
                [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_fields(Xmesh,Ymesh,nanowires,C,D,x_int,y_int,x_ext,y_ext);
                stop = toc;
                fprintf("\nIt took %.4f seconds to compute the fields.\n\n",stop)
                save(filename,'Etot','Htot','Einc','Hinc','Eref','Href','Escat','Hscat','X','Y','nanowires','C','D','x_int','y_int','x_ext','y_ext');
            end

            br = br + 1;
            waitbar(br/24,wb);

        end
    end
end

close(wb);

%% Compute relative error of scattered fields

clear; clc; close;

addpath('../../utils')

settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};

view = "far"; % close far far_approximation
choice = 3;

filename = sprintf('../../../../Results/nanowires/seperate_computation_%s_%s_N_100.mat',view,settings{choice});
load(filename);

predicted_field = Escat;

filename = sprintf('../../../../Results/nanowires/combined_computation_%s_%s_N_100.mat',view,settings{choice});
load(filename);

true_field = Escat;

relative = false;

error = compute_error(true_field,predicted_field,relative);
min(error,[],'all')
max(error,[],'all')


relative = true;
error = compute_error(true_field,predicted_field,relative);
min(error,[],'all')
max(error,[],'all')

filename = sprintf('../../../../Results/nanowires/seperate_computation_%s_%s_N_100.mat',view,settings{choice});
load(filename);

predicted_field = Hscat;

filename = sprintf('../../../../Results/nanowires/combined_computation_%s_%s_N_100.mat',view,settings{choice});
load(filename);

true_field = Hscat;

relative = false;

error = compute_error(true_field,predicted_field,relative);
min(error,[],'all')
max(error,[],'all')


relative = true;
error = compute_error(true_field,predicted_field,relative);
min(error,[],'all')
max(error,[],'all')

%% Setup test for two nanowires with some distance between them

clear; close all; clc;

choice = 5;
view = "far";
num_grid_points = 300;
N = 100;

distances = linspace(4*10^(-8),4*10^(-6),20);

dist = distances(end);

[nanowires, X, Y] = setup_settings(choice, N, num_grid_points, view, dist);

figure('Renderer', 'painters', 'Position', [400 400 1000 250]);
for j = 1:length(nanowires)
    nw = nanowires{j};
    plot(nw.x*10^6,nw.y*10^6,'.','Color',"#A2142F")
    hold on
    plot(nw.x_int*10^6,nw.y_int*10^6,'.','Color',"#0072BD")
    plot(nw.x_ext*10^6,nw.y_ext*10^6,'.','Color',"#77AC30")
end
axis equal
xlim([min(X)*10^6,max(X)*10^6])
ylim([-0.1,0.5])
grid on
xlabel('$x$ [$\mu$m]','FontSize',14);
ylabel('$y$ [$\mu$m]','FontSize',14);
t = linspace(-0.25,1.7,3);
plot(t,t*0,'k--')

absolute_E_error = 0 * distances;
relative_E_error = 0 * distances;
absolute_H_error = 0 * distances;
relative_H_error = 0 * distances;

wb = waitbar(0,'Computing fields...');
br = 0;
for j = 1:length(distances)
    dist = distances(j);
    [nanowires, X, Y] = setup_settings(choice, N, num_grid_points, view, dist);

    Escats = cell(2,1);
    Hscats = cell(2,1);
    
    for comp = 1:2
        
        if comp == 1
    
            % Run seperate computation
            tic;
            C = []; D = []; x_int = []; y_int = []; x_ext = []; y_ext = [];
            for k = 1:length(nanowires)
                temp = {nanowires{k}};
                [c,d,t_x_int,t_y_int,t_x_ext,t_y_ext] = forward(temp);
                C = [C c];
                D = [D d];
                x_int = [x_int; t_x_int];
                y_int = [y_int; t_y_int];
                x_ext = [x_ext; t_x_ext];
                y_ext = [y_ext; t_y_ext];
            end
            stop = toc;
            fprintf("\nIt took %.4f seconds to solve all the seperate linear systems.\n\n",stop)
    
        elseif comp == 2
    
            % Run combined computation
            tic;
            [C,D,x_int,y_int,x_ext,y_ext] = forward(nanowires);
            stop = toc;
            fprintf("\nIt took %.4f seconds to solve the combined linear system.\n\n",stop)
        end
    
        tic;
        [Etot, Htot, Escat, Hscat, Einc, Hinc, Eref, Href] = compute_far_fields(X,Y,C,x_int,y_int);
        stop = toc;
        fprintf("\nIt took %.4f seconds to compute the far fields.\n\n",stop)
    
        Escats{comp} = Escat;
        Hscats{comp} = Hscat;

        br = br + 1;
        waitbar(br/(2*length(distances)),wb);
    
    end
    
    predicted_field = Escats{1};
    true_field = Escats{2};
    
    relative = false;
    absolute_E_error(j) = max(compute_error(true_field,predicted_field,relative),[],'all');
    
    relative = true;
    relative_E_error(j) = max(compute_error(true_field,predicted_field,relative),[],'all');

    predicted_field = Hscats{1};
    true_field = Hscats{2};
    
    relative = false;
    absolute_H_error(j) = max(compute_error(true_field,predicted_field,relative),[],'all');
    
    relative = true;
    relative_H_error(j) = max(compute_error(true_field,predicted_field,relative),[],'all');

end

close(wb);

figure('Renderer', 'painters', 'Position', [400 400 1000 250]);
tiledlayout(1,3,'TileSpacing','compact');

nexttile;
plot(distances*10^6,absolute_E_error,'.-','LineWidth',1,'color',"#0072BD",'markersize',10)
grid on
ylabel('Absolute Error','fontsize', 14)
xlabel('Distance [$\mu$m]','fontsize', 14)

nexttile;
plot(distances*10^6,absolute_H_error,'o-','LineWidth',1,'color',"#D95319")
grid on
ylabel('Absolute Error','fontsize', 14)
xlabel('Distance [$\mu$m]','fontsize', 14)

nexttile;
plot(distances*10^6,relative_E_error*100,'.-','LineWidth',1,'color',"#0072BD",'DisplayName','\textbf{E}-field error','markersize',10)
hold on
plot(distances*10^6,relative_H_error*100,'o-','LineWidth',1,'color',"#D95319",'DisplayName','\textbf{H}-field error')
grid on
ylabel('Relative Error [\%]','fontsize', 14)
xlabel('Distance [$\mu$m]','fontsize', 14)
l = legend('NumColumns',2,'Box','off','fontsize',14);
l.Layout.Tile = "south";

destination = sprintf('../../../../Illustrations/nanowires/two_nanowires_with_distance_error.png');
exportgraphics(gcf,destination,'Resolution',300);

%%

clear; close all; clc;

addpath('../fields_2D')


settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
choice = 1;
filename = sprintf('../../../../Results/nanowires/combined_computation_close_%s_N_100.mat',settings{choice});
load(filename);

n = 1000;
[E_error, H_error] = transmission_error(C,D,x_int,y_int,x_ext,y_ext,nanowires,n);
for k = 1:length(nanowires)
    for hej = 1:14
        E_error(:,k) = smooth(E_error(:,k));
        H_error(:,k) = smooth(H_error(:,k));
    end
end
phi = linspace(0,2*pi,n);
num_nanowires = length(nanowires);

figure('Renderer', 'painters', 'Position', [400 400 800 250]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(phi,E_error,'LineWidth',1.5)
set(gca,'xtick',0:(pi/4):(2*pi), 'xticklabels',{'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}) 
ax=gca;
ax.XAxis.FontSize = 14;
grid on
xlim([0,2*pi])
xlabel('$\phi$','FontSize',14)
ylabel('Transmission Error','FontSize',14)
title('\textbf{E}-field','FontSize',14)

nexttile;
plot(phi,H_error,'LineWidth',1.5)
set(gca,'xtick',0:(pi/4):(2*pi), 'xticklabels',{'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}) 
ax=gca;
ax.XAxis.FontSize = 14;
grid on
xlim([0,2*pi])
title('\textbf{H}-field','FontSize',14)

xlabel('$\phi$','FontSize',14)

sgtitle(sprintf('\\textbf{Setup %d}',choice),'FontSize',16,'inteRpreter','latex')

destination = '../../../../Illustrations/nanowires/nanowire_transmission_error_one_nanowire.png';
exportgraphics(gcf,destination,'Resolution',300);





%%
computation = "seperate"; % seperate combined

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(2,3,'TileSpacing','compact');
for choice = 2:4
    
    filename = sprintf('../../../../Results/nanowires/%s_computation_close_%s_N_100.mat',computation,settings{choice});
    load(filename);
    
    n = 1000;
    [E_error, H_error] = transmission_error(C,D,x_int,y_int,x_ext,y_ext,nanowires,n);
    if strcmp(computation,"combined")
        for k = 1:length(nanowires)
            for hej = 1:14
                E_error(:,k) = smooth(E_error(:,k));
            end
        end
    end

    phi = linspace(0,2*pi,n);
    num_nanowires = length(nanowires);
    nexttile;
    if strcmp(computation,"combined")
        semilogy(phi,E_error,'LineWidth',1.5)
    else
        plot(phi,E_error,'LineWidth',1.5)
    end
    set(gca,'xtick',0:(pi/4):(2*pi), 'xticklabels',{'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}) 
    ax=gca;
    ax.XAxis.FontSize = 14;
    grid on
    xlim([0,2*pi])
    title(sprintf('\\textbf{Setup %d}',choice),'FontSize',14)
    xlabel('$\phi$','FontSize',14)
    if choice == 2
        ylabel({'\textbf{E}-field','Transmission Error'},'FontSize',14)
    end
    if strcmp(computation,"combined")
        ylim([10^(-10),10^(-4)])
    else
        ylim([0,1])
    end
end

for choice = 2:4
    
    filename = sprintf('../../../../Results/nanowires/%s_computation_close_%s_N_100.mat',computation,settings{choice});
    load(filename);
    
    n = 1000;
    [E_error, H_error] = transmission_error(C,D,x_int,y_int,x_ext,y_ext,nanowires,n);
    if strcmp(computation,"combined")
        for k = 1:length(nanowires)
            for hej = 1:14
                H_error(:,k) = smooth(H_error(:,k));
            end
        end
    end

    phi = linspace(0,2*pi,n);
    num_nanowires = length(nanowires);
    nexttile;
    if strcmp(computation,"combined")
        semilogy(phi,H_error,'LineWidth',1.5)
    else
        plot(phi,H_error,'LineWidth',1.5)
    end
    set(gca,'xtick',0:(pi/4):(2*pi), 'xticklabels',{'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}) 
    ax=gca;
    ax.XAxis.FontSize = 14;
    grid on
    xlim([0,2*pi])
    xlabel('$\phi$','FontSize',14)
    if choice == 2
        ylabel({'\textbf{H}-field','Transmission Error'},'FontSize',14)
    end
    if strcmp(computation,"combined")
        ylim([10^(-11),10^(-5)])
    else
        ylim([0,0.002])
    end
end

l = legend('Nanowire 1', 'Nanowire 2', 'Nanowire 3', 'Nanowire 4','box','off','numcolumns',4,'fontsize',14);
l.Layout.Tile = 'south';

if strcmp(computation,"combined")
    sgtitle('\textbf{Combined Computation}','fontsize',16,'interpreter','latex')
elseif strcmp(computation,"seperate")
    sgtitle('\textbf{Seperate Computation}','fontsize',16,'interpreter','latex')
end

destination = sprintf('../../../../Illustrations/nanowires/nanowire_transmission_error_%s_computation.png',computation);
exportgraphics(gcf,destination,'Resolution',300);
