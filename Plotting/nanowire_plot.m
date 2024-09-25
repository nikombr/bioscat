set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')

%% Plot results for one nanowire

clear; close all; clc;

settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
choice = 1;

views = {"far","close"}; % far_approximation

Etot_array = cell(2,1);
Htot_array = cell(2,1);
X_array = cell(2,1);
Y_array = cell(2,1);

for k = 1:2
    filename = sprintf('../Results/nanowires/seperate_computation_%s_%s_N_100.mat',views{k},settings{choice});
    load(filename);
    Etot_array{k} = Etot;
    Htot_array{k} = Htot;
    X_array{k}    = X;
    Y_array{k}    = Y;
end


climmaxE = cell(2,1);
climmaxH = cell(2,1);
for k = 1:2
    climmaxE{k} = -1;
    climmaxH{k} = -1;
    for i = 1:3
        climmaxE{k} = max(climmaxE{k},max(max(abs(Etot_array{k}(:,:,i)).^2)));
        climmaxH{k} = max(climmaxH{k},max(max(abs(Htot_array{k}(:,:,i)).^2)));
    end
end
climmaxE{1} = max(climmaxE{2},climmaxE{1});
climmaxE{2} = max(climmaxE{2},climmaxE{1});
climmaxH{1} = max(climmaxH{2},climmaxH{1});
climmaxH{2} = max(climmaxH{2},climmaxH{1});


figure('Renderer', 'painters', 'Position', [400 400 1000 450]);
tiledlayout(2,3,'TileSpacing','compact');
for k = 1:2
    Etot = Etot_array{k};
    Htot = Htot_array{k};
    X = X_array{k};
    Y = Y_array{k};

    nexttile;
    varnum = 3; % z
    if k == 1
        yshift = 3*10^(-2); % 3 cm
        ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
        plot_2D_field(X,Y,Etot,"E",climmaxE{k},varnum,yshift,ylabel)
    else
        plot_2D_field(X,Y,Etot,"E",climmaxE{k},varnum)
    end

    for j = 1:2
        nexttile;
        varnum = j; % z
        if k == 1
            yshift = 3*10^(-2); % 3 cm
            ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
            plot_2D_field(X,Y,Htot,"H",climmaxH{k},varnum,yshift,ylabel)
        else
            plot_2D_field(X,Y,Htot,"H",climmaxH{k},varnum)
        end
    end
end

sgtitle(sprintf('\\textbf{Setup %d}',choice),'interpreter','latex','fontsize',16)
destination = sprintf('../Illustrations/nanowires/fields_one_nanowire.png');
exportgraphics(gcf,destination,'Resolution',300);

%% Plot for two nanowires far from each other

clear; close all; clc;

settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
choice = 4;

views = {"far_approximation","close"};

computation = {"seperate", "combined"};

Etot_array = cell(2,2);
Htot_array = cell(2,2);
X_array = cell(2,2);
Y_array = cell(2,2);
reflectance_Ez = cell(2,1);
reflectance_Hx = cell(2,1);
%reflectance_Hy = cell(2,1);

% 2: {[0.051720804897082]} {[0.051723797644151]}
% 3: {[0.051720806196655]} {[0.052398680672409]}
% 4: {[0.050721077578794]} {[0.050054615075372]}

for j = 1:2
    for k = 1:2
        filename = sprintf('../Results/nanowires/%s_computation_%s_%s_N_100.mat',computation{j},views{k},settings{choice});
        load(filename);
        if k == 1
            reflectance_Ez{j} = mean((abs(Eref(:,:,3) + Escat(:,:,3))./abs(Einc(:,:,3))).^2,"all")
            % top = abs(Href + Hscat);
            % top = sqrt(top(:,:,1) + top(:,:,2));
            % bottom = abs(Hinc);
            % bottom = sqrt(bottom(:,:,1) + bottom(:,:,2));
            % reflectance_H{j} = mean((top./bottom).^2,"all")
            reflectance_Hx{j} = mean((abs(Href(:,:,1) + Hscat(:,:,1))./abs(Hinc(:,:,1))).^2,"all")
            %reflectance_Hy{j} = mean((abs(Href(:,:,2) + Hscat(:,:,2))./abs(Hinc(:,:,2))).^2,"all")
        end
        Etot_array{k,j} = Etot;
        Htot_array{k,j} = Htot;
        X_array{k,j}    = X;
        Y_array{k,j}    = Y;
    end
end

climmaxE = cell(2,1);
climmaxH = cell(2,1);
for j = 1:2
    for k = 1:2
        climmaxE{k} = -1;
        climmaxH{k} = -1;
        for i = 1:3
            climmaxE{k} = max(climmaxE{k},max(max(abs(Etot_array{k,j}(:,:,i)).^2)));
            climmaxH{k} = max(climmaxH{k},max(max(abs(Htot_array{k,j}(:,:,i)).^2)));
        end
    end
end
climmaxE{1} = max(climmaxE{2},climmaxE{1});
climmaxE{2} = max(climmaxE{2},climmaxE{1});
climmaxH{1} = max(climmaxH{2},climmaxH{1});
climmaxH{2} = max(climmaxH{2},climmaxH{1});

figure('Renderer', 'painters', 'Position', [400 400 1000 450]);
tiledlayout(2,3,'TileSpacing','compact');
varnum = 3; % z
for k = 1:2
    for j = 1:2
        Etot = Etot_array{k,j};
        Htot = Htot_array{k,j};
        X = X_array{k,j};
        Y = Y_array{k,j};
    
        nexttile;
        if k == 1
            yshift = 3*10^(-2); % 3 cm
            ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
            plot_2D_field(X,Y,Etot,"E",climmaxE{k},varnum,yshift,ylabel)
            if strcmp("seperate",computation{j})
                title('Seperate Computation','fontsize',14)
            elseif strcmp("combined",computation{j})
                title('Combined Computation','fontsize',14)
            end

        else
            plot_2D_field(X,Y,Etot,"E",climmaxE{k},varnum)
        end
    end
    
    nexttile;
    error = Etot_array{k,1}-Etot_array{k,2};
    if k == 1
        yshift = 3*10^(-2); % 3 cm
        ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
        plot_2D_field(X,Y,error,"E",-1,varnum,yshift,ylabel)
        title('Absolute Error','fontsize',14)
    else
        plot_2D_field(X,Y,error,"E",-1,varnum)
    end
    
end
sgtitle(sprintf('\\textbf{Setup %d}',choice),'interpreter','latex','fontsize',16)
destination = sprintf('../Illustrations/nanowires/Ez_field_%s.png',settings{choice});
exportgraphics(gcf,destination,'Resolution',300);


figure('Renderer', 'painters', 'Position', [400 400 1000 450]);
tiledlayout(2,3,'TileSpacing','compact');
varnum = 1; % x
for k = 1:2
    for j = 1:2
        Etot = Etot_array{k,j};
        Htot = Htot_array{k,j};
        X = X_array{k,j};
        Y = Y_array{k,j};
    
        nexttile;
        if k == 1
            yshift = 3*10^(-2); % 3 cm
            ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
            plot_2D_field(X,Y,Htot,"H",climmaxH{k},varnum,yshift,ylabel)
            if strcmp("seperate",computation{j})
                title('Seperate Computation','fontsize',14)
            elseif strcmp("combined",computation{j})
                title('Combined Computation','fontsize',14)
            end

        else
            plot_2D_field(X,Y,Htot,"H",climmaxH{k},varnum)
        end
    end
    
    nexttile;
    error = Htot_array{k,1}-Htot_array{k,2};
    if k == 1
        yshift = 3*10^(-2); % 3 cm
        ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
        plot_2D_field(X,Y,error,"H",-1,varnum,yshift,ylabel)
        title('Absolute Error','fontsize',14)
    else
        plot_2D_field(X,Y,error,"H",-1,varnum)
    end
    
end
sgtitle(sprintf('\\textbf{Setup %d}',choice),'interpreter','latex','fontsize',16)
destination = sprintf('../Illustrations/nanowires/Hx_field_%s.png',settings{choice});
exportgraphics(gcf,destination,'Resolution',300);

figure('Renderer', 'painters', 'Position', [400 400 1000 450]);
tiledlayout(2,3,'TileSpacing','compact');
varnum = 2; % x
for k = 1:2
    for j = 1:2
        Etot = Etot_array{k,j};
        Htot = Htot_array{k,j};
        X = X_array{k,j};
        Y = Y_array{k,j};
    
        nexttile;
        if k == 1
            yshift = 3*10^(-2); % 3 cm
            ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
            plot_2D_field(X,Y,Htot,"H",climmaxH{k},varnum,yshift,ylabel)
            if strcmp("seperate",computation{j})
                title('Seperate Computation','fontsize',14)
            elseif strcmp("combined",computation{j})
                title('Combined Computation','fontsize',14)
            end

        else
            plot_2D_field(X,Y,Htot,"H",climmaxH{k},varnum)
        end
    end
    
    nexttile;
    error = (Htot_array{k,1}-Htot_array{k,2});
    if k == 1
        yshift = 3*10^(-2); % 3 cm
        ylabel = '$y-3\textrm{ cm}$ [$\mu$m]';
        plot_2D_field(X,Y,error,"H",-1,varnum,yshift,ylabel)
        title('Absolute Error','fontsize',14)
    else
        plot_2D_field(X,Y,error,"H",-1,varnum)
    end
    
end
sgtitle(sprintf('\\textbf{Setup %d}',choice),'interpreter','latex','fontsize',16)
destination = sprintf('../Illustrations/nanowires/Hy_field_%s.png',settings{choice});
exportgraphics(gcf,destination,'Resolution',300);


%% Plot test points and auxiliary sources

clear; close all; clc;

colors = {"#0072BD","#D95319","#EDB120","#7E2F8E"};

settings = {'one_nanowire','two_far_nanowires','two_close_nanowires','multiple_nanowires'};
figure('Renderer', 'painters', 'Position', [400 400 1000 400]);
tiledlayout(2,2,'TileSpacing','compact');
for k = 1:4
    setting = settings{k};
    filename = sprintf('../Results/nanowires/combined_computation_far_%s_N_100.mat',setting);
    load(filename);

    nexttile;
    R = [];
    for j = 1:length(nanowires)
        nw = nanowires{j};
        c = [nw.xc nw.r]*10^6;
        r = nw.r*10^6;
        pos = [c-r 2*r 2*r];
        r = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', colors{j}, 'Edgecolor','none','FaceAlpha',0.2);
        hold on
        f = fill([10^(-8)],[10^(-8)],'k','Edgecolor', 'none','FaceColor', colors{j},'FaceAlpha',0.2);
        R = [R f];
    end
    for j = 1:length(nanowires)
        nw = nanowires{j};
        plot(nw.x*10^6,nw.y*10^6,'.','Color',"#A2142F")
        hold on
        plot(nw.x_int*10^6,nw.y_int*10^6,'.','Color',"#0072BD")
        plot(nw.x_ext*10^6,nw.y_ext*10^6,'.','Color',"#77AC30")
    end
    axis equal
    xlim([-0.25,1.7])
    ylim([-0.1,0.5])
    grid on
    xlabel('$x$ [$\mu$m]','FontSize',14);
    ylabel('$y$ [$\mu$m]','FontSize',14);
    t = linspace(-0.25,1.7,3);
    plot(t,t*0,'k--')
    title(sprintf('\\textbf{Setup %d}',k),'FontSize',14)
    
  


end

l = legend(R,{'Nanowire 1','Nanowire 2','Nanowire 3','Nanowire 4'},'NumColumns',4,'fontsize',14,'box','off');
l.Layout.Tile = 'south';





destination = '../Illustrations/nanowires/nanowire_setup_points.png';
exportgraphics(gcf,destination,'Resolution',300);
















