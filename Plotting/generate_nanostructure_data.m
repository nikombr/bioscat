set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')
clear; close all; clc;

filenames = {"Retinin2x2","demoleus2x2"};


figure('Renderer', 'painters', 'Position', [400 400 800 600]);
tiledlayout(2,2,'TileSpacing','compact');

for i = 1:2
    filename = filenames{i};

    data = readmatrix(sprintf("../Data/%s.xyz",filename),'FileType','text');
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    
    n = 1000;
    X = linspace(min(x),max(x),n);
    Y = linspace(min(y),max(y),1000);
    
    [valX,valY] = meshgrid(X,Y);
    
    val = griddata(x,y,z,valX,valY);

    if strcmp("Retinin2x2", filename)
        minimum_value = min(min(val));
    elseif strcmp("demoleus2x2", filename)
        val = val - min(min(val)) + minimum_value;
    end

    val = val + 10^(-8);
    
    nexttile
    imagesc(X,Y,val)
    colorbar
    title(filename)
    axis equal
    axis tight
    clim([0,max(max(val))])

    writematrix(valX, sprintf('../Data/nanostructures/3D/%s_x_%d.txt',filename,n));
    writematrix(valY, sprintf('../Data/nanostructures/3D/%s_y_%d.txt',filename,n));
    writematrix(val,  sprintf('../Data/nanostructures/3D/%s_f_%d.txt',filename,n));

    Y = val(550,:);
    X = X - mean(X);
    nexttile;
    plot(X,Y,'LineWidth',1.5)
    
    xlim([min(X),max(X)])
    %axis square
    grid on
    ylim([0,max(Y)])
    
    save(sprintf("../Data/%s_2D.mat",filename),"X", "Y")

    writematrix(X, sprintf('../Data/nanostructures/2D/%s_x_%d.txt',filename,n));
    writematrix(Y, sprintf('../Data/nanostructures/2D/%s_f_%d.txt',filename,n));
    


end







