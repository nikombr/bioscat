function preConditionedCrankNicholsonScheme(setup,protein_structure,total_grid_points,delta,covfunc,arg1,arg2,num_segments,data_quality)

% Set values
num = 10000;

if strcmp(protein_structure,"demoleus2x2")
    shift = 3*10^(-8);
elseif strcmp(protein_structure,"Retinin2x2")
    shift = 2*10^(-8);
end

% Choose functions to use
if strcmp(setup,"test")
    setupFunction = @testSetup;
    trueValueFunction = @testTrueValue;
    data_quality = "clean";
    num_segments = 1;

elseif strcmp(setup,"2D")
    setupFunction = @nanostructures2Dsetup;
    trueValueFunction = @nanostructures2DtrueValue;

elseif strcmp(setup,"3D")
    setupFunction = 0;
    trueValueFunction = 0;

end

%parpool;


% Get true values and arguments
[X, Y_true, true_val, args] = trueValueFunction(protein_structure, total_grid_points, num_segments, data_quality);

figure(1);
plot(X,Y_true,'k-','LineWidth',2)
hold on
xlim([min(X),max(X)])
grid on
destination = 'figure01.png';
exportgraphics(gcf,destination,'Resolution',300);
figure(4);
plot(X,Y_true,'k-','LineWidth',2)
hold on
xlim([min(X),max(X)])
grid on
destination = 'figure04.png';
exportgraphics(gcf,destination,'Resolution',300);

% Load the random curves
data = getRandomCurves(covfunc,arg1,arg2,X);

% Set initial curve
f = data(1,:)*0;
Y = f + shift;
n = 1;

% Get log-likelihood of initial curve
Lprev = setupFunction(X, Y, true_val, args{:});

alpha_array = zeros(2*num,1);
Lprev_array = zeros(2*num,1);
Lstar_array = zeros(2*num,1);
Lmean_array = zeros(2*num,1);
Y_array     = Y;
n_used = 0;
while n_used < num
    
    % Get random curve from file
    phi = data(n,:);

    % Compute possible new curve
    fstar = sqrt(1 - 2*delta)*f + sqrt(2*delta)*phi;

    % Shift to make sure that the curve is the right place
    Y = fstar + shift;

    % Plot values
    if n < 3000
        figure(1);
        destination = sprintf('%s/figure01.png',protein_structure);
    else
        figure(4);
        destination = sprintf('%s/figure04.png',protein_structure);
    end
    %plot(X,Y,'r-','LineWidth',0.3)
    %plot(X,phi,'y-','LineWidth',0.2)


    if mod(n,100) == 0
        plot(X,Y_true,'r-','LineWidth',2)
        plot(X,Y_mean,'m-','LineWidth',2)
        %plot(X,Y,'k-','LineWidth',2)
        %exportgraphics(gcf,destination,'Resolution',300);
    end
    
    % Compute log-likelihood
    Lstar = setupFunction(X, Y, true_val, args{:});
    
    Y_mean = mean(Y_array);

    % Compute temporary mean log-likelihood
    Lmean = setupFunction(X, Y_mean, true_val);
    
    % Compute probability of accepting curve
    alpha = min(1,exp(Lstar-Lprev));
    
    % Generate random number
    u = rand();

    % Determing whether or not to accept the curve
    if u < alpha
        f = fstar;
        Y_array = [Y_array; Y];
        
        Lprev = Lstar;
        plot(X,Y,'-','LineWidth',0.5, 'Color',[0.2 0.5 0.9 0.2])
        n_used = n_used + 1;
    else
        %plot(X,Y,'g-','LineWidth',0.3)
    end
        
    if mod(n,200) == 0
        figure(2);
        plot(alpha_array(1:n),'.-')
        title('$\alpha$')
        grid on
        destination = sprintf('%s/figure02.png',protein_structure);
        %exportgraphics(gcf,destination,'Resolution',300);

        figure(3);
        hold off
        plot(Lprev_array(1:n),'r.-')
        hold on
        plot(Lstar_array(1:n),'b.-')
        title('$L$')
        grid on
        destination = sprintf('%s/figure03.png',protein_structure);
        %exportgraphics(gcf,destination,'Resolution',300);

        figure(5);
        hold off
        semilogy(Lmean_array(1:n),'m.-')
        hold on
        title('$L$')
        grid on
        destination = sprintf('%s/figure05.png',protein_structure);
        %exportgraphics(gcf,destination,'Resolution',300);

        fprintf("The log-likelihood of the mean is %.4f\n\n",Lmean)
    end

    % Save values
    Lprev_array(n) = Lprev;
    Lstar_array(n) = Lstar;
    Lmean_array(n) = Lmean;
    alpha_array(n) = alpha;

    
    n = n + 1;



end

if strcmp(setup,"test")
    filename = sprintf('../../../Results/inverse_%s/%s_%s_%s_%s_total_x_grid_points_%d.mat',setup,covfunc,arg1,arg2,protein_structure,total_grid_points);
else
    filename = sprintf('../../../Results/inverse_%s/%s/%s_%s_%s_%s_total_x_grid_points_%d_num_segments_%d.mat',setup,data_quality,covfunc,arg1,arg2,protein_structure,total_grid_points,num_segments);
end

save(filename,"Lprev_array","Lstar_array","alpha_array","Y_array","X","Y_true")