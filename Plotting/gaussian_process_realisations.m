set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')


%% Curves

clear; close all; clc;

ells = [1 2 4];
taus = {"0.25", "0.5", "1"};

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(3,3,'TileSpacing','compact');
for k = 1:length(ells)
    for j = 1:length(taus)
        ell = ells(k);
        tau = taus{j};
        filename = sprintf("../Data/gaussian_process_realisations/curve_squared_exponential_tau_%s_ell_%d.mat",tau,ell);
        load(filename)

        [m,n] = size(data);

        nexttile;
        for i = 1:5
            plot(x,data(i,:))
   
            hold on
        end
        ylim([-2,2])
        grid on
        title(sprintf("$\\tau=%s$, $\\ell=%d$",tau,ell),'FontSize',16)
    end


end
save_folder = "../Illustrations/gaussian_process_realisations/";
exportgraphics(gcf,sprintf('%scurve_squared_exponential_n_%d.png',save_folder,n),'Resolution',300);

%% Planes



clear; close all; clc;

ells = [1 2 4];
taus = {"0.25", "0.5", "1"};

minclim = 10000;
maxclim = -10000;

for k = 1:length(ells)
    for j = 1:length(taus)
        ell = ells(k);
        tau = taus{j};
        filename = sprintf("../Data/gaussian_process_realisations/plane_squared_exponential_tau_%s_ell_%d.mat",tau,ell);
        load(filename)
        minclim = min(minclim, min(min(data)));
        maxclim = max(maxclim, max(max(data)));
    end

end



for k = 1:length(ells)
    figure('Renderer', 'painters', 'Position', [400 400 1000 600]);
    tiledlayout(3,5,'TileSpacing','compact');
    for j = 1:length(taus)
        ell = ells(k);
        tau = taus{j};
        filename = sprintf("../Data/gaussian_process_realisations/plane_squared_exponential_tau_%s_ell_%d.mat",tau,ell);
        load(filename)

        [m,n] = size(data);
        n = length(x);

        for i = 1:5
            nexttile;
            imagesc(x,y,reshape(data(i,:),[n,n]))
   
            hold on
            grid on
            if i == 3
                title(sprintf("$\\tau=%s$",tau),'FontSize',16)
            end
            axis square
            clim([minclim, maxclim])
            
            
        end
        
    end
    c = colorbar;
    c.Layout.Tile = 'east';
     sgtitle(sprintf("$\\ell=%d$",ell),'FontSize',16,'interpreter','latex')
     save_folder = "../Illustrations/gaussian_process_realisations/";
     exportgraphics(gcf,sprintf('%splane_squared_exponential_ell_%d_n_%d.png',save_folder,ell,n),'Resolution',300);

end
