function preConditionedCrankNicholsonScheme()
parpool;
dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';

num_segments = 10;
total_x_grid_points = 1000;
protein_structure = 'demoleus2x2';
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

load(sprintf('%sData/gaussian_process_realisations/curve_matern_p_1_ell_2.mat',dir),'data','x');


protein_structure = 'backward';

filename = sprintf("%sData/backward_2D.mat",dir);

minx = 4.8828e-09;
maxx = 2.0068e-06;

X = x/max(x)*(maxx-minx) + minx;

f = data(1,:);
%f = f - min(f);
f = f*10^(-8) + 3*10^(-8);

Y = f;

save(filename,"Y","X");

plot(X,Y,'m-','LineWidth',1.5)

[RE_last] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas);

delta = 0.3;

num = 30;
alpha_array = [];
Lprev_array = [];

for n = 2:num
    
    
    phi = data(n,:);
    %phi = phi - min(phi);
    phi = phi*10^(-8) + 3*10^(-8);

    fstar = sqrt(1-2*delta)*f+sqrt(2*delta)*phi;
    Y = fstar;
    figure(1);
    plot(X,Y,'r-','LineWidth',0.5)
    plot(X,phi,'y-','LineWidth',0.5)
    
    save(filename,"Y","X");
    
    tic;
    [RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas);
    stop = toc;
    fprintf("\nIt took %.4f seconds to compute the reflectance.\n\n",stop)

    Lstar = -sum(abs(RE_true-RE).^2,'all')/2;
    Lprev = -sum(abs(RE_true-RE_last).^2,'all')/2;
    Lprev_array = [Lprev_array Lprev];

    alpha = min(1,exp(Lstar-Lprev));
    alpha_array = [alpha_array alpha];

    u = rand();
    figure(1);
    if u < alpha
        f = fstar;
        RE_last = RE;
        plot(X,Y,'b-','LineWidth',0.5)
    else
        plot(X,Y,'g-','LineWidth',0.5)
    end
    figure(2);
    plot(alpha_array,'.-')
    grid on
    title('$\alpha$')

    figure(3);
    plot(Lprev_array,'.-')
    grid on
    title('$L$')



end

alpha_array