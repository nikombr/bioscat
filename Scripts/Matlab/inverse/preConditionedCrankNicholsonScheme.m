function preConditionedCrankNicholsonScheme()
%parpool;
dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
dir = '/zhome/00/b/147112/bioscat/';

num_segments = 20;
total_x_grid_points = 1000;
protein_structure = 'Retinin2x2';  % Retinin2x2 demoleus2x2
data_quality = 'clean';

filename = sprintf('%sData/reflectance_2D/%s/%s_total_x_grid_points_%d_num_segments_%d.mat',dir,data_quality,protein_structure,total_x_grid_points,num_segments);
load(filename,'RE','coord_obs','lambdas','betas')


filename = sprintf("%sData/%s_2D.mat",dir,protein_structure);
load(filename,'X','Y')

X_true = X;
Y_true = Y;

figure(1);
plot(X_true,Y_true,'k-','LineWidth',2)
hold on
xlim([min(X),max(X)])
grid on
destination = 'figure01.png';
exportgraphics(gcf,destination,'Resolution',300);
figure(4);
plot(X_true,Y_true,'k-','LineWidth',2)
hold on
xlim([min(X),max(X)])
grid on
destination = 'figure04.png';
exportgraphics(gcf,destination,'Resolution',300);
RE_true = RE;


load(sprintf('%sData/gaussian_process_realisations/curve_matern_p_1_ell_0.2.mat',dir),'data','x');

protein_structure = 'backward';

filename = sprintf("%sData/backward_2D.mat",dir);

minx = 4.8828e-09;
maxx = 2.0068e-06;
minx = min(X_true);
maxx = max(X_true);

X = x/max(x)*(maxx-minx) + minx;

n = 1;
f = 0;
while n < 2 %min(f) < 10^(-9)
    f = data(n,:);
    f = f*10^(-8);%+ 2*10^(-8);
    n = n + 1;
end
%f = f - min(f);

%f = interp1(X_true,Y_true,X);
Y = f;
Y_last = Y;
%save(filename,"Y","X");

[RE_last] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y);
size(RE_last)
delta = 0.3;

num = 10000;
alpha_array = [];
Lprev_array = [];
Lstar_array = [];
Lmean_array = [];
F_array     = [f];

n_used = 1;

while n < num
    
    
    phi = data(n,:);
    %phi = phi - min(phi);
    phi = phi*10^(-8); %+ 10^(-7);
    minphi = min(phi);
    fstar = sqrt(1-2*delta)*f + sqrt(2*delta)*phi;
    
 
    %if min(fstar) > 10^(-9)
        
        
        
        %fstar = sqrt(1-2*delta)*f+sqrt(2*delta)*phi;
        Y = fstar + 2*10^(-8);
        if n_used < 3000
            figure(1);
            destination = 'figure01.png';
        else
            figure(4);
            destination = 'figure04.png';
        end
        plot(X,Y,'r-','LineWidth',0.3)
        plot(X,phi,'y-','LineWidth',0.2)
        plot(X,mean(F_array) + 2*10^(-8),'m-','LineWidth',2)
        
        if mod(n_used,20) == 0
            plot(X_true,Y_true,'k-','LineWidth',2)
            %plot(X,Y,'k-','LineWidth',2)
            exportgraphics(gcf,destination,'Resolution',300);
        end
    
        
        %save(filename,"Y","X");
        
        %start = tic;
        %[RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y);
        %stop = toc(start);
        %fprintf("\nIt took %.4f seconds to compute the reflectance.\n\n",stop)
        %Lstar = -sqrt(sum(abs(RE_true-RE).^2,'all'))/2;
        %size(Y_true)
        %size(Y)
        %size(Y_last)
        min(Y_true);
        max(Y_true);
        Lstar = -sum(abs((Y_true-Y)*10^8).^2,'all');
        Lprev = -sum(abs((Y_true-Y_last)*10^8).^2,'all');
        Lmean = -sum(abs((Y_true-(mean(F_array) + 2*10^(-8)))*10^8).^2,'all');
        %Lprev = -sqrt(sum(abs(RE_true-RE_last).^2,'all'))/2;
        Lprev_array = [Lprev_array Lprev];
        Lstar_array = [Lstar_array Lstar];
        Lmean_array = [Lmean_array Lmean];
    
        alpha = min(1,exp(Lstar-Lprev));
        alpha_array = [alpha_array alpha];
    
        u = rand();
        if n_used < 3000
            figure(1);
            destination = 'figure01.png';
        else
            figure(4);
            destination = 'figure04.png';
        end
        if n_used == 1000
            "deleted"
            F_array = [f];
        end
        if u < alpha
            f = fstar;
            F_array = [F_array; f];
            
            %size(F_array)
            RE_last = RE;
            Y_last = Y;
            plot(X,Y,'b-','LineWidth',0.5)
        else
            plot(X,Y,'g-','LineWidth',0.3)
        end
        
        %exportgraphics(gcf,destination,'Resolution',300);
        if mod(n_used,20) == 0
            figure(2);
            plot(alpha_array,'.-')
            title('$\alpha$')
            grid on
            destination = 'figure02.png';
            exportgraphics(gcf,destination,'Resolution',300);
    
            figure(3);
            hold off
            plot(Lprev_array,'r.-')
            hold on
            plot(Lstar_array,'b.-')
            title('$L$')
            grid on
            destination = 'figure03.png';
            exportgraphics(gcf,destination,'Resolution',300);

            figure(5);
            hold off
            semilogy(Lmean_array,'m.-')
            hold on
            title('$L$')
            grid on
            destination = 'figure05.png';
            exportgraphics(gcf,destination,'Resolution',300);
        end
        n_used = n_used + 1;
    %end
    n = n + 1;



end
