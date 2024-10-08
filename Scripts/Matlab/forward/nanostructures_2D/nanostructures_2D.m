
clear; close all; clc;

addpath('../fields_2D/')

protein_structure    = "Retinin2x2"; % Retinin2x2 demoleus2x2
num_segments_choices = [1 5 10 20]; % 1 5 10 20
total_x_grid_points  = 1000;

views = {"far","close","far_approximation"};

% Illustrate nanostructure
load(sprintf("../../../../Data/%s_2D.mat",protein_structure))
figure('Renderer', 'painters', 'Position', [400 400 800 250]);
plot(X,Y,'LineWidth',1.5)
xlim([min(X),max(X)])
grid on
hold on

% Do precomputations
x = linspace(min(X),max(X),total_x_grid_points);
d = x(2)-x(1);
alpha = 2*d;

moveY = 1.05*max(Y);

coord = struct;

%wb = waitbar(0,'Computing fields...');
%br = 0;
for num_segments = num_segments_choices
    load(sprintf("../../../../Data/%s_2D.mat",protein_structure))
    n = round(total_x_grid_points/num_segments);
    m = num_segments * n;

    x = linspace(min(X),max(X),m + 1);
    y = interp1(X,Y,x);

    segments = cell(num_segments,1);

    % Setup segments
    for k = 1:num_segments
    
        segx = x((k-1)*n+1:k*n+1);
        segy = y((k-1)*n+1:k*n+1);
        segments{k} = setup_nanostructures(segx,segy,alpha);

    end

    % Solve linear system for each segment
    tic;
    for k = 1:length(segments)
        segments{k} = forward(segments{k});
    end
    stop = toc;
    n = 300;
    fprintf("\nIt took %.4f seconds to solve all the linear systems.\n\n",stop)

    for view_choice = 1:3
        view = views{view_choice};

        % Form grid in which the fields are computed in
        Y = linspace(0,21*10^(-7),n);
        X = linspace(-0.5*10^(-7),20.5*10^(-7),n);
        if strcmp("far",view) || strcmp("far_approximation",view)
            Y = Y + 3*10^(-2);
        else
            Y = Y + moveY;
        end
        [Xmesh,Ymesh] = meshgrid(X,Y);
        coord.x = Xmesh;
        coord.y = Ymesh;
        
        % Compute the fields
        tic;
        if strcmp(view,"far_approximation")
            [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(coord, segments);
        else
            [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_fields(coord, segments);
        end
        stop = toc;
        fprintf("\nIt took %.4f seconds to compute the fields.\n\n",stop)
        
        % Save results
        filename = sprintf('../../../../Results/nanostructures_2D/%s_%s_num_segments_%d.mat',protein_structure,view,num_segments);
        save(filename,'Etot','Htot','Einc','Hinc','Eref','Href','Escat','Hscat','X','Y','segments');

        %br = br + 1;
        %waitbar(br/(3*length(num_segments_choices)),wb);
    end
end
%close(wb);
%%
figure;
imagesc(x,y,abs(Etot(:,:,3)).^2)
set(gca,'YDir','normal')
colorbar
hold on
for j = 1:length(segments{1}.x)
    %plot(segments{1}.x{j},segments{1}.y{j},'r.','LineWidth',1.5,'Markersize',10)
end
axis equal
axis tight

%% Illustrate the segments

clear; close all; clc;

view = "close";
num_segments = 5;
protein_structure = "demoleus2x2"; % Retinin2x2 demoleus2x2

filename = sprintf('../../../../Results/nanostructures_2D/%s_%s_num_segments_%d.mat',protein_structure,view,num_segments);
load(filename);


figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
for k = 1:num_segments
    
    segx = x((k-1)*n+1:k*n+1);
    segy = y((k-1)*n+1:k*n+1);
    
    segments{k} = setup_nanostructures(segx,segy,alpha);
    for j =1:length(segments{k}.x)
        plot(segments{k}.x{j},segments{k}.y{j},'r.','LineWidth',1.5,'Markersize',10)
        hold on
    end
    
    plot(segments{k}.x_int,segments{k}.y_int,'b.','LineWidth',1.5,'Markersize',10)
    
    plot(segments{k}.x_ext,segments{k}.y_ext,'go','LineWidth',1.5,'Markersize',5)
    plot(segments{k}.n_x,segments{k}.n_y,'k.','LineWidth',1.5,'Markersize',10)
end
xlim([-0.5*10^(-7),20.5*10^(-7)])
ylim([-0.5*10^(-7),5*10^(-7)])
grid on
axis equal

%% Compute relative error of scattered fields

clear; clc; close;
protein_structure = "Retinin2x2"; % Retinin2x2 demoleus2x2

addpath('../../utils')

num_segments = 20;

view = "far"; % close far far_approximation
choice = 3;

filename = sprintf('../../../../Results/nanostructures_2D/%s_%s_num_segments_%d.mat',protein_structure,view,num_segments);
load(filename);

predicted_field = Escat;

filename = sprintf('../../../../Results/nanostructures_2D/%s_%s_num_segments_1.mat',protein_structure,view);
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

filename = sprintf('../../../../Results/nanostructures_2D/%s_%s_num_segments_%d.mat',protein_structure,view,num_segments);
load(filename);

predicted_field = Hscat;

filename = sprintf('../../../../Results/nanostructures_2D/%s_%s_num_segments_1.mat',protein_structure,view);
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

%% Compute error for varying number of segments and track CPU time as well

clear; close all; clc;

%num_grid_points = 300;
%N = 100;

% Set values
total_x_grid_points  = 1000;
num_segments_array = 2:20;
protein_structure = "Retinin2x2"; % Retinin2x2 demoleus2x2
show_waitbar = false;

% Load data
load(sprintf("../../../../Data/%s_2D.mat",protein_structure))

% Form grid in which the fields are computed in
n = 100;
Ygrid = linspace(0,21*10^(-7),n) + 3*10^(-2);
Xgrid = linspace(-0.5*10^(-7),20.5*10^(-7),n);
[Xgrid, Ygrid] = meshgrid(Xgrid, Ygrid);
            
% Do precomputations
x = linspace(min(X),max(X),total_x_grid_points);
d = x(2)-x(1);
alpha = 2*d;

% Compute true solution
num_segments = 1;
n = round(total_x_grid_points/num_segments);
m = num_segments * n;

x = linspace(min(X),max(X),m + 1);
y = interp1(X,Y,x);

segments = cell(num_segments,1);

% Setup segments
for k = 1:num_segments
    segx = x((k-1)*n+1:k*n+1);
    segy = y((k-1)*n+1:k*n+1);
    segments{k} = setup_nanostructures(segx,segy,alpha);
end

CPUtime = zeros(length(num_segments_array)+1,1);

% Solve linear system for each segment
tic;
for k = 1:length(segments)
    segments{k} = forward(segments{k});
end
stop = toc;
CPUtime(1) = stop;
fprintf("\nIt took %.4f seconds to solve all the linear system.\n\n",stop)

tic;
[Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(Xgrid,Ygrid,segments,show_waitbar);
stop = toc;
fprintf("\nIt took %.4f seconds to compute the far fields.\n\n",stop)

% Save as true field
true_E_field = Escat;
true_H_field = Hscat;

absolute_E_error = 0 * num_segments_array;
relative_E_error = 0 * num_segments_array;
absolute_H_error = 0 * num_segments_array;
relative_H_error = 0 * num_segments_array;


wb = waitbar(0,'Computing fields and errors...');
br = 0;
for j = 1:length(num_segments_array)
    num_segments = num_segments_array(j);

    % Compute predicted solution
    n = round(total_x_grid_points/num_segments);
    m = num_segments * n;
    
    x = linspace(min(X),max(X),m + 1);
    y = interp1(X,Y,x);
    
    segments = cell(num_segments,1);
    
    % Setup segments
    for k = 1:num_segments
        segx = x((k-1)*n+1:k*n+1);
        segy = y((k-1)*n+1:k*n+1);
        segments{k} = setup_nanostructures(segx,segy,alpha);
    end
    
    % Solve linear system for each segment
    tic;
    for k = 1:length(segments)
        segments{k} = forward(segments{k});
    end
    stop = toc;
    CPUtime(j+1) = stop;
    fprintf("\nIt took %.4f seconds to solve all the linear system.\n\n",stop)
    
    tic;
    [Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(Xgrid,Ygrid,segments,show_waitbar);
    stop = toc;
    fprintf("\nIt took %.4f seconds to compute the far fields.\n\n",stop)

    predicted_E_field = Escat;
    predicted_H_field = Hscat;        
    
    relative = false;
    absolute_E_error(j) = max(compute_error(true_E_field,predicted_E_field,relative),[],'all');
    
    relative = true;
    relative_E_error(j) = max(compute_error(true_E_field,predicted_E_field,relative),[],'all');
    
    relative = false;
    absolute_H_error(j) = max(compute_error(true_H_field,predicted_H_field,relative),[],'all');
    
    relative = true;
    relative_H_error(j) = max(compute_error(true_H_field,predicted_H_field,relative),[],'all');

    br = br + 1;
    waitbar(br/length(num_segments_array),wb);

end

close(wb);

figure('Renderer', 'painters', 'Position', [400 400 1000 250]);
tiledlayout(1,3,'TileSpacing','compact');

nexttile;
plot(num_segments_array,absolute_E_error,'.-','LineWidth',1,'color',"#0072BD",'markersize',10)
grid on
ylabel('Absolute Error','fontsize', 14)
xlabel('Number of Segments','fontsize', 14)

nexttile;
plot(num_segments_array,absolute_H_error,'o-','LineWidth',1,'color',"#D95319")
grid on
ylabel('Absolute Error','fontsize', 14)
xlabel('Number of Segments','fontsize', 14)

nexttile;
plot(num_segments_array,relative_E_error*100,'.-','LineWidth',1,'color',"#0072BD",'DisplayName','\textbf{E}-field error','markersize',10)
hold on
plot(num_segments_array,relative_H_error*100,'o-','LineWidth',1,'color',"#D95319",'DisplayName','\textbf{H}-field error')
grid on
ylabel('Relative Error [\%]','fontsize', 14)
xlabel('Number of Segments','fontsize', 14)
l = legend('NumColumns',2,'Box','off','fontsize',14);
l.Layout.Tile = "south";

destination = sprintf('../../../../Illustrations/nanostructures_2D/%s_error_variying_number_of_segments.png',protein_structure);
exportgraphics(gcf,destination,'Resolution',300);


figure('Renderer', 'painters', 'Position', [400 400 400 250]);

plot(1:num_segments_array(end),CPUtime,'.-','LineWidth',1.5,'color',"#0072BD",'markersize',15)
grid on
ylabel('CPU time [s]','fontsize', 14)
xlabel('Number of Segments','fontsize', 14)

destination = sprintf('../../../../Illustrations/nanostructures_2D/%s_CPUtime.png',protein_structure);
exportgraphics(gcf,destination,'Resolution',300);


%% Far-field patterns

clear; clc; close all; 

addpath('../../utils')
addpath('../fields_2D')

% Set values
total_x_grid_points  = 1000;
num_segments_array = 2:20;
protein_structure = "demoleus2x2"; % Retinin2x2 demoleus2x2
show_waitbar = false;

% Load data
load(sprintf("../../../../Data/%s_2D.mat",protein_structure))
            
% Do precomputations
x = linspace(min(X), max(X), total_x_grid_points);
d = x(2) - x(1);
alpha = 2*d;
dim = 2; % 2D computation
num = 500; % Number of points to compute far-field pattern
show_waitbar = true;
r = 3*10^(-2); % Radius of where to compute far field pattern
num_segments_array = [1 5 10 15 20];
PE = cell(length(num_segments_array),1);
PH = cell(length(num_segments_array),1);

for j = 1:length(num_segments_array)
    
    % Compute solution
    num_segments = num_segments_array(j);
    n = round(total_x_grid_points/num_segments);
    m = num_segments * n;
    
    x = linspace(min(X),max(X),m + 1);
    y = interp1(X,Y,x);
    
    segments = cell(num_segments,1);
    
    % Setup segments
    for k = 1:num_segments
        segx = x((k-1)*n+1:k*n+1);
        segy = y((k-1)*n+1:k*n+1);
        segments{k} = setup_nanostructures(segx,segy,alpha);
    end
    
    % Solve linear system for each segment
    tic;
    for k = 1:length(segments)
        segments{k} = forward(segments{k});
    end
    stop = toc;
    fprintf("\nIt took %.4f seconds to solve all the linear system.\n\n",stop)
    
    tic;
    [phi, PE{j}, PH{j}] = compute_far_field_pattern(dim,num,r,segments,show_waitbar);
    %[Etot, Htot, Einc, Hinc, Eref, Href, Escat, Hscat] = compute_far_fields(Xgrid,Ygrid,);
    stop = toc;
    fprintf("\nIt took %.4f seconds to compute the far-field patterns.\n\n",stop)

end

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(2,2,'TileSpacing','compact');

for j = 2:5
    nexttile
    polarplot(phi,PE{1},'-','linewidth',1,'DisplayName','Combined Solution')
    hold on 
    polarplot(phi,PE{j},'linewidth',1.5,'DisplayName','Sum of Seperate Solutions')
    thetalim([0 180])
    if strcmp(protein_structure,'Retinin2x2')
        rlim([0,0.025])
    elseif strcmp(protein_structure,'demoleus2x2')
        rlim([0,0.03])
    end
    title(sprintf("%d Segments",num_segments_array(j)),'fontsize',14)
    thetaticks([0 30 60 90 120 150 180])
    thetaticklabels({'$0$','$\pi/6$','$\pi/3$','$\pi/2$','$2\pi/3$','$5\pi/6$','$\pi$'})
    %text(3*pi/2,0.005,'Intensity')
end

l = legend('numcolumns',2,'FontSize',14,'box','off');
l.Layout.Tile = 'south';

destination = sprintf('../../../../Illustrations/nanostructures_2D/%s_far_field_patterns.png',protein_structure);
exportgraphics(gcf,destination,'Resolution',300);

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(2,2,'TileSpacing','compact');
for j = 2:5
    nexttile
    polarplot(phi,abs(PE{1}-PE{j})/PE{1}*100,'linewidth',1.5)
    thetalim([0 180])
    
    title(sprintf("%d Segments",num_segments_array(j)),'fontsize',14)
    if strcmp(protein_structure,'Retinin2x2')
        rlim([0,10])
        rticks(0:3:9)
        rticklabels({'0\,\%','3\,\%','6\,\%','9\,\%'})
    elseif strcmp(protein_structure,'demoleus2x2')
        rlim([0,16])
        rticks(0:4:16)
        rticklabels({'0\,\%','4\,\%','8\,\%','12\,\%','16\,\%'})
    end
    
    thetaticks([0 30 60 90 120 150 180])
    thetaticklabels({'$0$','$\pi/6$','$\pi/3$','$\pi/2$','$2\pi/3$','$5\pi/6$','$\pi$'})
    %text(3*pi/2,0.005,'Intensity')
end


destination = sprintf('../../../../Illustrations/nanostructures_2D/%s_far_field_patterns_relative_error.png',protein_structure);
exportgraphics(gcf,destination,'Resolution',300);