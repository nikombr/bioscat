clear; close all; clc;

setup = "test";
total_grid_points = 1000;
protein_structure = 'demoleus2x2';  % Retinin2x2 demoleus2x2
covfunc = "matern";

deltas = [0.001 0.01 0.1 0.2 0.3];
arg1 = {"1","2","4"};
arg2 = {"0.1","0.2","0.4","0.8"};

delta = deltas(1);
num = 7000;


figure('Renderer', 'painters', 'Position', [400 400 1000 1000]);
tiledlayout(4,3,'TileSpacing','compact');
for j = 1:length(arg2)
    for i = 1:length(arg1)
        fprintf("(i, j) = (%d, %d)\n",i,j)
        filename = sprintf('../../../Results/inverse_%s/%s_%s_%s_%s_total_x_grid_points_%d_delta_%.3f_num_%d.mat',setup,protein_structure,covfunc,arg1{i},arg2{j},total_grid_points,delta,num);
        
            load(filename)
            Y_mean = mean(Y_array(1:3500,:));
            Lmean = testSetup(X, Y_mean, Y_true)
            nexttile;
            plot(X,Y_array(1:3500,:),'-','linewidth',0.1,'Color',[0 0.4470 0.7410 0.03])
            hold on
            plot(X,Y_mean,'b-','linewidth',2)
            plot(X,Y_true,'k-','linewidth',1)
            grid on
            xlabel('$x$','fontsize',14)
            
            xlim([min(X),max(X)])
            ylim([0,5*10^(-8)])
            if j == 1
                title({sprintf("$p=%d$",str2num(arg1{i})),sprintf("$L=%.4f$",Lmean)},'fontsize',15)
            else
                title(sprintf("$L=%.4f$",Lmean),'fontsize',15)
            end
            if i == 1
                ylabel({sprintf('$\\ell=%.1f$',str2num(arg2{j})),'$f(x)$'},'fontsize',14)
            else
                ylabel('$f(x)$','fontsize',14)
            end
        
    end
end
"hej"
destination = sprintf('%s/inverse_test_plot_first_half_delta_%.3f.png',protein_structure,delta);
exportgraphics(gcf,destination,'Resolution',300);

figure('Renderer', 'painters', 'Position', [400 400 1000 1000]);
tiledlayout(4,3,'TileSpacing','compact');
for j = 1:length(arg2)
    for i = 1:length(arg1)
        fprintf("(i, j) = (%d, %d)\n",i,j)
        filename = sprintf('../../../Results/inverse_%s/%s_%s_%s_%s_total_x_grid_points_%d_delta_%.3f_num_%d.mat',setup,protein_structure,covfunc,arg1{i},arg2{j},total_grid_points,delta,num);
        
            load(filename)
            Y_mean = mean(Y_array(3501:7000,:));
            Lmean = testSetup(X, Y_mean, Y_true)
            nexttile;
            plot(X,Y_array(3501:7000,:),'-','linewidth',0.1,'Color',[0 0.4470 0.7410 0.03])
            hold on
            plot(X,Y_mean,'b-','linewidth',2)
            plot(X,Y_true,'k-','linewidth',1)
            grid on
            xlabel('$x$','fontsize',14)
            ylabel('$f(x)$','fontsize',14)
            xlim([min(X),max(X)])
            ylim([0,5*10^(-8)])
            if j == 1
                title({sprintf("$p=%d$",str2num(arg1{i})),sprintf("$L=%.4f$",Lmean)},'fontsize',15)
            else
                title(sprintf("$L=%.4f$",Lmean),'fontsize',15)
            end
            if i == 1
                ylabel({sprintf('$\\ell=%.1f$',str2num(arg2{j})),'$f(x)$'},'fontsize',14)
            else
                ylabel('$f(x)$','fontsize',14)
            end
   
        
    end
end
"hej"
destination = sprintf('%s/inverse_test_plot_second_half_delta_%.3f.png',protein_structure,delta);
exportgraphics(gcf,destination,'Resolution',300);