set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')

%% Plot some results fast for quick visualisation
clear; close all; clc; 

filename = "demoleus2x2"; % Retinin2x2 demoleus2x2

climmaxE = -1;
climmaxHx = -1;
climmaxHy = -1;

segment_numbers = [1 5 10 20];

for num_segments = segment_numbers

    load(sprintf("../Results/nanostructures_2D/%s_num_segments_%d.mat",filename,num_segments))
    Etot = Etot(20:end,:,:);
    Htot = Htot(20:end,:,:);
    for k = 1:3
        climmaxE = max(climmaxE,max(max(abs(Etot(:,:,k)).^2)));
    end
    climmaxHx = max(climmaxHx,max(max(abs(Htot(:,:,1)).^2)));
    climmaxHy = max(climmaxHy,max(max(abs(Htot(:,:,2)).^2)));

end


figure('Renderer', 'painters', 'Position', [400 400 1400 900]);
tiledlayout(3,4,'TileSpacing','compact');

for num_segments = segment_numbers

    load(sprintf("../Results/nanostructures_2D/%s_num_segments_%d.mat",filename,num_segments))
    varnum = 3;
    Etot = Etot(20:end,:,:);
    Htot = Htot(20:end,:,:);
    X = X(20:end);
    Y = Y(20:end);
    nexttile;
    plot_2D_field(X,Y,Etot,"E",climmaxE,varnum)

    if num_segments == 1
        title(sprintf('%d segment',num_segments),'FontSize',14)
    else
        title(sprintf('%d segments',num_segments),'FontSize',14)
    end
end

for num_segments = segment_numbers

    load(sprintf("../Results/nanostructures_2D/%s_num_segments_%d.mat",filename,num_segments))
    varnum = 1;
    Etot = Etot(20:end,:,:);
    Htot = Htot(20:end,:,:);
    X = X(20:end);
    Y = Y(20:end);

    nexttile;
    plot_2D_field(X,Y,Htot,"H",climmaxHx,varnum)

end

for num_segments = segment_numbers

    load(sprintf("../Results/nanostructures_2D/%s_num_segments_%d.mat",filename,num_segments))
    varnum = 2;
    Etot = Etot(20:end,:,:);
    Htot = Htot(20:end,:,:);
    X = X(20:end);
    Y = Y(20:end);
    nexttile;
    plot_2D_field(X,Y,Htot,"H",climmaxHy,varnum)

end

destination = sprintf('../Illustrations/nanostructures_2D/%s_test.png',filename);
exportgraphics(gcf,destination,'Resolution',300);

%% Plot results with one segment

clear; close all; clc; 

filename = 'Retinin2x2'; % Retinin2x2 demoleus2x2
views = {'far', 'close'};


climmaxE = -1;
climmaxHx = -1;
climmaxHy = -1;

for j = 1:2
    
    view = views{j};
    load(sprintf("../Results/nanostructures_2D/%s_%s_num_segments_1.mat",filename,view))
    Etot = Etot(50:end,:,:);
    Htot = Htot(50:end,:,:);
    for k = 1:3
        climmaxE = max(climmaxE,max(max(abs(Etot(:,:,k)).^2)));
    end
    climmaxHx = max(climmaxHx,max(max(abs(Htot(:,:,1)).^2)));
    climmaxHy = max(climmaxHy,max(max(abs(Htot(:,:,2)).^2)));
end

    

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(2,3,'TileSpacing','compact');

for k = 1:2
    view = views{k};
    load(sprintf("../Results/nanostructures_2D/%s_%s_num_segments_1.mat",filename,view))

    

    if k == 2
        Etot = Etot(50:end,:,:);
        Htot = Htot(50:end,:,:);
        X = X(50:end);
        Y = Y(50:end);
    end

    if k == 1
        yshift = 3*10^(-2); % 3 cm
        ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
    else
        yshift = 0;
        ylabel = '$y$ [$\mu$m]';
    end
    
    varnum = 3;
    nexttile;
    plot_2D_field(X,Y,Etot,"E",climmaxE,varnum,yshift,ylabel)

    varnum = 1;
    nexttile;
    plot_2D_field(X,Y,Htot,"H",climmaxHx,varnum,yshift,ylabel)

    varnum = 2;
    nexttile;
    plot_2D_field(X,Y,Htot,"H",climmaxHy,varnum,yshift,ylabel)



end

destination = sprintf('../Illustrations/nanostructures_2D/%s_1_segment.png',filename);
exportgraphics(gcf,destination,'Resolution',300);

%% Plot error 

clear; close all; clc; 

filename = 'demoleus2x2'; % Retinin2x2 demoleus2x2
views = {'far', 'close'};

num_segments = 5;

climmaxE = -1;
climmaxHx = -1;
climmaxHy = -1;

for j = 1:2
    
    view = views{j};
    load(sprintf("../Results/nanostructures_2D/%s_%s_num_segments_1.mat",filename,view))
    Etot = Etot(50:end,:,:);
    Htot = Htot(50:end,:,:);
    for k = 1:3
        climmaxE = max(climmaxE,max(max(abs(Etot(:,:,k)).^2)));
    end
    climmaxHx = max(climmaxHx,max(max(abs(Htot(:,:,1)).^2)));
    climmaxHy = max(climmaxHy,max(max(abs(Htot(:,:,2)).^2)));
end

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
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

    load(sprintf("../Results/nanostructures_2D/%s_%s_num_segments_1.mat",filename,view))

    if k == 2
        Etot_1 = Etot(50:end,:,:);
        Htot_1 = Htot(50:end,:,:);
        X = X(50:end);
        Y = Y(50:end);
    else
        Etot_1 = Etot;
    end
    
    varnum = 3;

    nexttile;
    plot_2D_field(X,Y,Etot_1,"E",climmaxE,varnum,yshift,ylabel)

    load(sprintf("../Results/nanostructures_2D/%s_%s_num_segments_%d.mat",filename,view,num_segments))

    if k == 2
        Etot = Etot(50:end,:,:);
        Htot = Htot(50:end,:,:);
        X = X(50:end);
        Y = Y(50:end);
    end
 
    nexttile;
    plot_2D_field(X,Y,Etot,"E",climmaxE,varnum,yshift,ylabel)

    nexttile;
    plot_2D_field(X,Y,Etot-Etot_1,"E",-1,varnum,yshift,ylabel)



end

destination = sprintf('../Illustrations/nanostructures_2D/%s_Ez_%d_segments.png',filename,num_segments);
exportgraphics(gcf,destination,'Resolution',300);