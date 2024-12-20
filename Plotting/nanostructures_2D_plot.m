
%% Plot results with one segment

clear; close all; clc; 

protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2
views = {'far', 'near'};
scenario = 1;

if scenario == 1
    beta = 0;
elseif scenario == 2
    beta = 90;
end
lambda=325; % nm
total_grid_points = 300;
climmaxE = -1;
climmaxH = -1;


for j = 1:2
    
    view = views{j};
    %%load(sprintf("../Results/nanostructures_2D/%s_scenario_%d_%s_num_segments_1.mat",filename,scenario,view))
    filename = sprintf("../Results/forward/%s/%s/fields_beta_%d_lambda_%d_num_segments_1_total_grid_points_%d.mat",protein_structure,view,beta,lambda,total_grid_points);
    load(filename)
    E_tot = E_ref + E_scat;
    H_tot = H_ref + H_scat;
    for k = 1:3
        climmaxE = max(climmaxE,max(max(abs(E_tot(:,:,k)).^2)));
        climmaxH = max(climmaxH,max(max(abs(H_tot(:,:,k)).^2)));
    end
    
end

if scenario == 1
    fieldtype = {"E","H","H"};
    varnums = [3 1 2];
elseif scenario == 2
    fieldtype = {"E","E","H"};
    varnums = [1 2 3];
end

figure('Renderer', 'painters', 'Position', [400 400 1000 550]);
tiledlayout(2,3,'TileSpacing','compact');

for k = 1:2
    view = views{k};
    %load(sprintf("../Results/nanostructures_2D/%s_scenario_%d_%s_num_segments_1.mat",filename,scenario,view))
    filename = sprintf("../Results/forward/%s/%s/fields_beta_%d_lambda_%d_num_segments_1_total_grid_points_%d.mat",protein_structure,view,beta,lambda,total_grid_points);
    load(filename)
    E_tot = E_ref + E_scat;
    H_tot = H_ref + H_scat;
    if k == 1
        yshift = 3*10^(-2); % 3 cm
        ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
    else
        yshift = 0;
        ylabel = '$y$ [$\mu$m]';
    end
    for j = 1:3
        if strcmp(fieldtype{j}, "E")
            field = E_tot;
            climmax = climmaxE;
        elseif strcmp(fieldtype{j}, "H")
            field = H_tot;
            climmax = climmaxH;
        end
        varnum = varnums(j);
        nexttile;
        plot_2D_field(x,y,field,fieldtype{j},climmax,varnum,yshift,ylabel)
    end

end
annotation_of_type(scenario,protein_structure)

destination = sprintf('../Illustrations/nanostructures_2D/%s_scenario_%d_1_segment.png',protein_structure,scenario);
exportgraphics(gcf,destination,'Resolution',300);

%% Plot error 

clear; close all; clc; 

filename = 'Retinin2x2'; % Retinin2x2 demoleus2x2
views = {'far', 'close'};

for num_segments = [5 10 20]

for scenario = 1:2

vars = {'x', 'y', 'z'};

if scenario == 1
    fieldtype = {"E","H","H"};
    varnums = [3 1 2];
elseif scenario == 2
    fieldtype = {"E","E","H"};
    varnums = [1 2 3];
end

climmaxE = -1;
climmaxH = -1;

for j = 1:2
    
    view = views{j};
    load(sprintf("../Results/nanostructures_2D/%s_scenario_%d_%s_num_segments_1.mat",filename,scenario,view))
    for k = 1:3
        climmaxE = max(climmaxE,max(max(abs(Etot(:,:,k)).^2)));
        climmaxH = max(climmaxH,max(max(abs(Htot(:,:,k)).^2)));
    end

    load(sprintf("../Results/nanostructures_2D/%s_scenario_%d_%s_num_segments_%d.mat",filename,scenario,view,num_segments))
    for k = 1:3
        climmaxE = max(climmaxE,max(max(abs(Etot(:,:,k)).^2)));
        climmaxH = max(climmaxH,max(max(abs(Htot(:,:,k)).^2)));
    end
end

for j = 1:3

    varnum = varnums(j);

    figure('Renderer', 'painters', 'Position', [400 400 1000 550]);
    tiledlayout(2,3,'TileSpacing','compact');
    
    for k = 1:2
        view = views{k};
        
        if k == 1
            yshift = 3*10^(-2); % 3 cm
            ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
        else
            yshift = 0;
            ylabel = '$y$ [$\mu$m]';
        end
    
        load(sprintf("../Results/nanostructures_2D/%s_scenario_%d_%s_num_segments_1.mat",filename,scenario,view))
        
        if strcmp(fieldtype{j},"E")
            field_1 = Etot;
            climmax = climmaxE;
        elseif strcmp(fieldtype{j},"H")
            field_1 = Htot;
            climmax = climmaxH;
        end
        
    
        nexttile;
        plot_2D_field(X,Y, field_1,fieldtype{j},climmax,varnum,yshift,ylabel)
        if k == 1
            title('1 Segment','fontsize',14)
        end
    
        load(sprintf("../Results/nanostructures_2D/%s_scenario_%d_%s_num_segments_%d.mat",filename,scenario,view,num_segments))

        if strcmp(fieldtype{j},"E")
            field = Etot;
        elseif strcmp(fieldtype{j},"H")
            field = Htot;
        end
        
       
     
        nexttile;
        plot_2D_field(X,Y,field, fieldtype{j},climmax,varnum,yshift,ylabel)
        
        if k == 1
            title(sprintf('%d Segments',num_segments),'fontsize',14)
        end
    
        nexttile;
        plot_2D_field(X,Y,field - field_1, fieldtype{j},-1,varnum,yshift,ylabel)
        if k == 1
            title('Absolute Error','fontsize',14)
        end
    
    
    end

    annotation_of_type(scenario,filename)
    
    destination = sprintf('../Illustrations/nanostructures_2D/%s_scenario_%d_%s%s_%d_segments.png',filename,scenario,fieldtype{j},vars{varnum},num_segments);
    exportgraphics(gcf,destination,'Resolution',300);

end

close all;
end
end


%% Illustrate reflectance


clear; close all; clc;

data_quality = 'noisy'; % clean noisy


protein_structure    = "demoleus2x2"; % Retinin2x2 demoleus2x2
num_segments = 20; % 1 5 10 20
total_x_grid_points  = 1000;

for num_segments = [1 5 10 20]
% Load clean data (without noise)
filename = sprintf('../Data/reflectance_2D/%s/%s_total_x_grid_points_%d_num_segments_%d.mat',data_quality,protein_structure,total_x_grid_points,num_segments);
load(filename)

if num_segments == 1
    RE_true = RE;
else
    fprintf('The relative error from using %d segments instead of 1 is %.2f %%\n\n',num_segments,max(abs(RE_true-RE)./abs(RE_true),[],'all')*100);
end

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(2,3,'TileSpacing','compact');
for k = 1:6
    nexttile;
    for j = 1:7
        plot(betas,RE(j,:,k),'DisplayName',sprintf('$\\lambda=%.1f$ nm',lambdas(j)*10^(9)),'LineWidth',1.5)
        hold on
    end
    ylim([3.5,5.3])
    grid on
    exp_x = floor(log10(abs(coord_obs.x(k))));
    exp_y = floor(log10(abs(coord_obs.y(k))));
    xlabel('$\beta$ [radians]','fontsize',14);
    ylabel('Reflectance','fontsize',14);
    title(sprintf('$x=%.1f\\cdot10^{%d}$ m, $y=%.1f\\cdot10^{%d}$ m',coord_obs.x(k)*10^(-exp_x),exp_x,coord_obs.y(k)*10^(-exp_y),exp_y))
end
l = legend('box','off','numcolumns',7,'fontsize',14);
l.ItemTokenSize = [20 20];
l.Layout.Tile = "south";

sgtitle(sprintf('\\textbf{Number of Segments: %d}',num_segments),'fontsize',14,'interpreter','latex')

annotation_of_type(-1,protein_structure,'protein_type')

destination = sprintf('../Illustrations/nanostructures_2D/reflectance/%s_%s_total_x_grid_points_%d_num_segments_%d.png',protein_structure,data_quality,total_x_grid_points,num_segments);
exportgraphics(gcf,destination,'Resolution',300);

end