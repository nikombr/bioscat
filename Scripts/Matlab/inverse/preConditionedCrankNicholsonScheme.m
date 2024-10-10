function preConditionedCrankNicholsonScheme()
parpool;
dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
dir = '/zhome/00/b/147112/bioscat/';

num_segments = 10;
total_x_grid_points = 1000;
protein_structure = 'Retinin2x2';  % Retinin2x2 demoleus2x2
data_quality = 'noisy';

filename = sprintf('%sData/reflectance_2D/%s/%s_total_x_grid_points_%d_num_segments_%d.mat',dir,data_quality,protein_structure,total_x_grid_points,num_segments);
load(filename)


filename = sprintf("%sData/%s_2D.mat",dir,protein_structure);
load(filename)

figure(1);
plot(X,Y,'k-','LineWidth',2)
hold on
xlim([min(X),max(X)])
grid on
RE_true = RE;

X_hmm = X;
Y_hmm = Y;

load(sprintf('%sData/gaussian_process_realisations/curve_matern_p_4_ell_0.5.mat',dir),'data','x');


protein_structure = 'backward';

filename = sprintf("%sData/backward_2D.mat",dir);

minx = 4.8828e-09;
maxx = 2.0068e-06;

X = x/max(x)*(maxx-minx) + minx;

f = data(2,:);
%f = f - min(f);
f = f*10^(-8) + 2*10^(-8);
f = interp1(X_hmm,Y_hmm,X);
Y = f;

save(filename,"Y","X");

plot(X,Y,'m-','LineWidth',1.5)

[RE_last] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas);

delta = 0.2;

num = 2000;
alpha_array = [];
Lprev_array = [];

for n = 3:num
    
    
    phi = data(n,:);
    %phi = phi - min(phi);
    phi = phi*10^(-8) + 2*10^(-8);

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
    destination = 'figure01.png';
    if mod(n,10)
        %plot(X,Y,'k-','LineWidth',2)
        exportgraphics(gcf,destination,'Resolution',300);
    end

    
    save(filename,"Y","X");
    
    tic;
    [RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas);
    stop = toc;
    fprintf("\nIt took %.4f seconds to compute the reflectance.\n\n",stop)
    RE_true
    RE_last
    Lstar = -sum(abs(RE_true-RE).^2,'all')/2;
    Lprev = -sum(abs(RE_true-RE_last).^2,'all')/2
    Lprev_array = [Lprev_array Lprev];

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
        plot(alpha_array)
        destination = 'figure02.png';
        exportgraphics(gcf,destination,'Resolution',300);

        figure(3);
        plot(Lprev_array)
        destination = 'figure03.png';
        exportgraphics(gcf,destination,'Resolution',300);
    end
    
    



end
