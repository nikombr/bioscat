function preConditionedCrankNicholsonScheme(setup,num,protein_structure,total_grid_points,delta,covfunc,arg1,arg2,num_segments,data_quality)

% Set values
%num = 7000;

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
Y_true = 10^(-8)*sin(X*10^7) + shift;
true_val = compute_reflectance("backward", args{:}, X, Y_true);
figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
plot(X,Y_true,'k-','LineWidth',2)
hold on
xlim([min(X),max(X)])
grid on
destination = 'figure01.png';
exportgraphics(gcf,destination,'Resolution',300);
figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
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

maxnum = 20000;

alpha_array = zeros(maxnum,1);
Lprev_array = zeros(maxnum,1);
Lstar_array = zeros(maxnum,1);
Lmean_array = zeros(maxnum,1);
Y_array     = zeros(num, total_grid_points);
Y_array(1,:) = Y;
Y_mean = Y;
n_used = 1;

print_frequency = 10;

while n_used < num && n < maxnum
    
    % Get random curve from file
    phi = data(n,:);

    % Compute possible new curve
    fstar = sqrt(1 - 2*delta)*f + sqrt(2*delta)*phi;

    % Shift to make sure that the curve is the right place
    Y = fstar + shift;

    % Plot values
    if mod(n,print_frequency) == 0
        if n < 3000
            figure(1);
            destination = sprintf('%s/figure01.png',protein_structure);
        else
            figure(2);
            destination = sprintf('%s/figure02.png',protein_structure);
        end
        plot(X,Y_true,'r-','LineWidth',2)
        plot(X,Y_mean,'b-','LineWidth',2)
        %plot(X,Y,'k-','LineWidth',2)
        exportgraphics(gcf,destination,'Resolution',300);
        plot(X,Y_mean,'m-','LineWidth',2)
    end
    
    % Compute log-likelihood
    Lstar = setupFunction(X, Y, true_val, args{:});
    
    % Compute probability of accepting curve
    alpha = min(1,exp(Lstar-Lprev));
    
    % Generate random number
    u = rand();

    % Determing whether or not to accept the curve
    if u < alpha
        f = fstar;
        
        Y_array(n_used+1,:) = Y;
        
        Lprev = Lstar;
        plot(X,Y,'-','LineWidth',1, 'Color',[0.2 0.5 0.9 0.2])
        n_used = n_used + 1;
    else
        %plot(X,Y,'g-','LineWidth',0.3)
    end
        
    if mod(n,print_frequency) == 0
        if n_used >2000
            Y_mean = mean(Y_array(2000:n_used,:));
        else
            Y_mean = mean(Y_array(1:n_used,:));
        end

        % Compute temporary mean log-likelihood
        Lmean = setupFunction(X, Y_mean, true_val, args{:})
        Lmean = testSetup(X, Y_mean, Y_true);

       %Lmean = 0;

        figure(4);
        plot(alpha_array(1:n),'.-')
        title('$\alpha$')
        grid on
        destination = sprintf('%s/figure04.png',protein_structure);
        exportgraphics(gcf,destination,'Resolution',300);

        figure(3);
        hold off
        plot(Lprev_array(1:n),'r.-')
        hold on
        plot(Lstar_array(1:n),'b.-')
        title('$L$')
        grid on
        destination = sprintf('%s/figure03.png',protein_structure);
        exportgraphics(gcf,destination,'Resolution',300);

        figure(5);
        semilogy(n,abs(Lmean),'m.-','markersize',20)
        hold on
        title('$L$')
        grid on
        destination = sprintf('%s/figure05.png',protein_structure);
        exportgraphics(gcf,destination,'Resolution',300);

        fprintf("The log-likelihood of the mean is %.4f and %d have been used.\n\n",Lmean,n_used)
    end

    % Save values
    Lprev_array(n) = Lprev;
    Lstar_array(n) = Lstar;
    %Lmean_array(n) = Lmean;
    alpha_array(n) = alpha;

    
    n = n + 1;



end

% Scale so we do not store too much data
Y_array_new = zeros(num,500);
X_new = linspace(min(X),max(X),500);
for i = 1:num
    Y_array_new(i,:) = interp1(X,Y_array(i,:),X_new);
end
Y_true = interp1(X,Y_true,X_new);
X = X_new;
Y_array = Y_array_new;

if strcmp(setup,"test")
    filename = sprintf('../../../Results/inverse_%s/%s_%s_%s_%s_total_x_grid_points_%d_delta_%.3f_num_%d.mat',setup,protein_structure,covfunc,arg1,arg2,total_grid_points,delta,num);
else
    filename = sprintf('../../../Results/inverse_%s/%s/%s_%s_%s_%s_total_x_grid_points_%d_num_segments_%d_delta_%.3f_num_%d.mat',setup,protein_structure,data_quality,covfunc,arg1,arg2,total_grid_points,num_segments,delta,num);
end

save(filename,"Lprev_array","Lstar_array","alpha_array","Y_array","X","Y_true","n_used","n")