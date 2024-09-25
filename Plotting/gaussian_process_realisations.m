set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex')


%% Curves

clear; close all; clc;

covfunc = "squared_exponential"; % "matern"

ells = [1 2 4];
taus = {"0.25", "0.5", "1"};
ps = [1 2 4];

maxylim = -1;

for k = 1:3
    for j = 1:3
        ell = ells(k);
        tau = taus{j};
        p = ps(j);
        if strcmp(covfunc, "squared_exponential")
            filename = sprintf("../Data/gaussian_process_realisations/curve_%s_tau_%s_ell_%d.mat",covfunc,tau,ell);
        elseif strcmp(covfunc, "matern")
            filename = sprintf("../Data/gaussian_process_realisations/curve_%s_p_%d_ell_%d.mat",covfunc,p,ell);
        end
        load(filename)
        maxylim = max(maxylim, abs(max(max(data(1:5,:)))));
    end

end

maxylim = ceil(maxylim);

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(3,3,'TileSpacing','compact');
for k = 1:3
    for j = 1:3
        ell = ells(k);
        tau = taus{j};
        p = ps(j);
        if strcmp(covfunc, "squared_exponential")
            filename = sprintf("../Data/gaussian_process_realisations/curve_%s_tau_%s_ell_%d.mat",covfunc,tau,ell);
        elseif strcmp(covfunc, "matern")
            filename = sprintf("../Data/gaussian_process_realisations/curve_%s_p_%d_ell_%d.mat",covfunc,p,ell);
        end
        load(filename)

        [m,n] = size(data);

        nexttile;
        for i = 1:5
            plot(x,data(i,:),'LineWidth',1)
   
            hold on
        end
        ylim([-maxylim,maxylim])
        grid on
        if strcmp(covfunc, "squared_exponential")
             title(sprintf("$\\tau=%s$, $\\ell=%d$",tau,ell),'FontSize',16)
        elseif strcmp(covfunc, "matern")
             title(sprintf("$p=%d$, $\\ell=%d$",p,ell),'FontSize',16)
        end
    end


end
save_folder = "../Illustrations/gaussian_process_realisations/";
exportgraphics(gcf,sprintf('%scurve_%s_n_%d.png',save_folder,covfunc,n),'Resolution',300);

%% Planes



clear; close all; clc;

covfunc = "matern"; % "matern"  "squared_exponential"
ps = [1 2 4];

ells = [1 2 4];
taus = {"0.25", "0.5", "1"};

minclim = 10000;
maxclim = -10000;

for k = 1:length(ells)
    for j = 1:length(taus)
        ell = ells(k);
        tau = taus{j};
        p = ps(j);
        if strcmp(covfunc, "squared_exponential")
            filename = sprintf("../Data/gaussian_process_realisations/plane_%s_tau_%s_ell_%d.mat",covfunc,tau,ell);
        elseif strcmp(covfunc, "matern")
            filename = sprintf("../Data/gaussian_process_realisations/plane_%s_p_%d_ell_%d.mat",covfunc,p,ell);
        end
        load(filename)
        minclim = min(minclim, min(min(data(1:5,:))));
        maxclim = max(maxclim, max(max(data(1:5,:))));
    end

end



for k = 1:length(ells)
    figure('Renderer', 'painters', 'Position', [400 400 1000 600]);
    tiledlayout(3,5,'TileSpacing','compact');
    for j = 1:length(taus)
        ell = ells(k);
        tau = taus{j};
        p = ps(j);
        if strcmp(covfunc, "squared_exponential")
            filename = sprintf("../Data/gaussian_process_realisations/plane_%s_tau_%s_ell_%d.mat",covfunc,tau,ell);
        elseif strcmp(covfunc, "matern")
            filename = sprintf("../Data/gaussian_process_realisations/plane_%s_p_%d_ell_%d.mat",covfunc,p,ell);
        end
        load(filename)

        [m,n] = size(data);
        n = length(x);

        for i = 1:5
            nexttile;
            imagesc(x,y,reshape(data(i,:),[n,n]))
   
            hold on
            
            if i == 3
                if strcmp(covfunc, "squared_exponential")
                    title(sprintf("$\\tau=%s$",tau),'FontSize',16)
                elseif strcmp(covfunc, "matern")
                    title(sprintf("$p=%d$",p),'FontSize',16)
                end
                
            end
            axis square
            clim([minclim, maxclim])
            
            
        end
        
    end
    c = colorbar;
    c.Layout.Tile = 'east';
     sgtitle(sprintf("$\\ell=%d$",ell),'FontSize',16,'interpreter','latex')
     save_folder = "../Illustrations/gaussian_process_realisations/";
     exportgraphics(gcf,sprintf('%splane_%s_ell_%d_n_%d.png',save_folder,covfunc,ell,n),'Resolution',300);

end
