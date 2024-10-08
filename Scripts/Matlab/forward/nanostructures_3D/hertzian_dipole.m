clear; clc; % close all; 

addpath("../fields_3D")

x_ext = 0.5*10^(-7);
%x_ext = 0;
y_ext = 0;
z_ext = 0;

grid_num = 300;
Zvalues = [-2*10^(-7) -10^(-7) 0 10^(-7) 2*10^(-7)];
Ztitles = {"-2\cdot10^{-7}","-10^{-7}","0","10^{-7}","2\cdot10^{-7}"};


figure('Renderer', 'painters', 'Position', [400 400 800 1200]);
tiledlayout(5,3,'TileSpacing','compact');

for j = 1:5
    Y = linspace(-2,2,grid_num)*10^(-7);
    X = linspace(-2,2,grid_num)*10^(-7);
    [Xmesh,Ymesh] = meshgrid(X,Y);
    
    % Get size of field to compute
    [n,m] = size(Xmesh);
    Xmesh = reshape(Xmesh,[],1);
    Ymesh = reshape(Ymesh,[],1);
    M = n*m;
    Zmesh = Xmesh*0 + Zvalues(j);
    
    E = E_tot_inside_matrix(Xmesh,Ymesh,Zmesh,x_ext,y_ext,z_ext);
    E = reshape(E,n,m,3);
    
    for k = 1:3
        nexttile;
        imagesc(X,Y,abs(E(:,:,k)).^2)
        colorbar;
        axis equal
        axis tight
        if k == 2
            title(sprintf("$z = %s$ m",Ztitles{j}),'FontSize',14)
        end
    end
end

Xvalues = [-2*10^(-7) -10^(-7) 0 10^(-7) 2*10^(-7)];
Xtitles = {"-2\cdot10^{-7}","-10^{-7}","0","10^{-7}","2\cdot10^{-7}"};


destination = '../../../../Illustrations/nanostructures_3D/hertzian_dipole_z.png';
exportgraphics(gcf,destination,'Resolution',300);


figure('Renderer', 'painters', 'Position', [400 400 800 1200]);
tiledlayout(5,3,'TileSpacing','compact');

for j = 1:5
    Y = linspace(-2,2,grid_num)*10^(-7);
    Z = linspace(-2,2,grid_num)*10^(-7);
    [Ymesh,Zmesh] = meshgrid(Y,Z);
    
    % Get size of field to compute
    [n,m] = size(Ymesh);
    Ymesh = reshape(Ymesh,[],1);
    Zmesh = reshape(Zmesh,[],1);
    M = n*m;
    Xmesh = Ymesh*0 + Xvalues(j);
    
    E = E_tot_inside_matrix(Xmesh,Ymesh,Zmesh,x_ext,y_ext,z_ext);
    E = reshape(E,n,m,3);
    
    for k = 1:3
        nexttile;
        imagesc(X,Y,abs(E(:,:,k)).^2)
        colorbar;
        axis equal
        axis tight
        if k == 2
            title(sprintf("$x = %s$ m",Xtitles{j}),'FontSize',14)
        end
    end
end

destination = '../../../../Illustrations/nanostructures_3D/hertzian_dipole_x.png';
exportgraphics(gcf,destination,'Resolution',300);