function preConditionedCrankNicholsonScheme()
%parpool;
dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
%dir = '/zhome/00/b/147112/bioscat/';

num_segments = 10;
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
RE_true = RE;


load(sprintf('%sData/gaussian_process_realisations/curve_matern_p_4_ell_0.8.mat',dir),'data','x');

protein_structure = 'backward';

filename = sprintf("%sData/backward_2D.mat",dir);

minx = 4.8828e-09;
maxx = 2.0068e-06;
minx = min(X_true);
maxx = max(X_true);

X = x/max(x)*(maxx-minx) + minx;

f = data(2,:);
%f = f - min(f);
f = f*10^(-8) + 2*10^(-8);
%f = interp1(X_true,Y_true,X);
Y = f;

%save(filename,"Y","X");

%plot(X,Y,'m-','LineWidth',1.5)

[RE_last] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y);

delta = 0.001;

num = 2000;
alpha_array = [];
Lprev_array = [];
Lstar_array = [];

for n = 3:num
    
    n
    phi = data(n,:);
    %phi = phi - min(phi);
    phi = phi*10^(-8) + 2*10^(-8);
    minphi = min(phi)
    if all(phi > 10^(-8))
    
    

    fstar = sqrt(1-2*delta)*f+sqrt(2*delta)*phi;
    Y = fstar;
    if n < 1000
        figure(1);
        destination = 'figure01.png';
    else
        figure(4);
        destination = 'figure04.png';
    end
    plot(X,Y,'r-','LineWidth',0.5)
    plot(X,phi,'y-','LineWidth',0.5)
    %destination = 'figure01.png';
    if mod(n,10)
        %plot(X,Y,'k-','LineWidth',2)
        exportgraphics(gcf,destination,'Resolution',300);
    end

    
    %save(filename,"Y","X");
    
    tic;
    [RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y);
    stop = toc;
    fprintf("\nIt took %.4f seconds to compute the reflectance.\n\n",stop)
    RE_last
    Lstar = -sum(abs(RE_true-RE).^2,'all')*4;
    Lprev = -sum(abs(RE_true-RE_last).^2,'all')*4;
    Lprev_array = [Lprev_array Lprev];
    Lstar_array = [Lstar_array Lstar];

    alpha = min(1,exp(Lstar-Lprev));
    alpha_array = [alpha_array alpha];

    u = rand();
    if n < 1000
        figure(1);
        destination = 'figure01.png';
    else
        figure(4);
        destination = 'figure04.png';
    end
    if u < alpha
        f = fstar;
        RE_last = RE;
        plot(X,Y,'b-','LineWidth',0.5)
    else
        plot(X,Y,'g-','LineWidth',0.5)
    end
    
    %exportgraphics(gcf,destination,'Resolution',300);
    if mod(n,10)
        figure(2);
        plot(alpha_array,'.-')
        title('$\alpha$')
        %destination = 'figure02.png';
        %exportgraphics(gcf,destination,'Resolution',300);

        figure(3);
        hold off
        plot(Lprev_array,'.-')
        hold on
        plot(Lstar_array,'.-')
        title('$L$')
        %destination = 'figure03.png';
        %exportgraphics(gcf,destination,'Resolution',300);
    end
    end
    



end
