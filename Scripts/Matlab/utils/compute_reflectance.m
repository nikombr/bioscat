function [RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas, X, Y)
% beta: angle of displacement of polarisation
% lambda0: wavelength of incident plane wave in free space
% coord_obs: where we want to observe the reflectance


far_field_approximation = false; %

if strcmp(protein_structure,'backward')
    % If we are doing inverse/backward computations, we do not save
    % segments to a file
    segments = setup_segments(X,Y,num_segments,total_x_grid_points);
else
    dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';
    dir = '/zhome/00/b/147112/bioscat/';
    segment_filename = sprintf('%sData/segments_2D/%s_total_x_grid_points_%d_num_segments_%d.mat',dir,protein_structure,total_x_grid_points,num_segments);
    if ~isfile(segment_filename)
        % Setup segments if they are not yet saved
        setup_settings(protein_structure, total_x_grid_points, num_segments);
    else
        % Get how the segments are planned out
        load(segment_filename,'segments');
    end
end



% Loop sizes
n1 = length(lambdas);
n2 = 2;
n3 = length(segments);

% Allocate for segments
compute_segments = cell(n1*n2*n3,1);

% Move variables to processes
lambdas        = parallel.pool.Constant(lambdas);
%local_segments = parallel.pool.Constant(local_segments);
segments       = parallel.pool.Constant(segments);
betas = parallel.pool.Constant(betas);
coord_obs = parallel.pool.Constant(coord_obs);

% Do forward for all segments needed
%tic;
%ticBytes(gcp);
parfor idx = 1:(n1*n2*n3)

    % Suppresses singular matrix warnings
    warning('off', 'MATLAB:singularMatrix');  
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:rankDeficientMatrix')

    % Get true indices
    k = mod(floor((idx-1) / (n2 * n3)), n1) + 1;
    scenario = mod(floor((idx-1) / n3), n2) + 1;
    j = mod(idx-1, n3) + 1;

    % Get value
    lambda = lambdas.Value(k);

    % Do forward
    compute_segments{idx} = forward(segments.Value{j},scenario,lambda,false);
  
end
%tocBytes(gcp)
%stop = toc;

%fprintf('It took %.4f seconds to comupute the segments.\n',stop)
%tic;
%ticBytes(gcp);
parfor k = 1:n1
    M = length(coord_obs.Value.x);
    RE_local = zeros(length(betas.Value),M);

    far_field_approximation = false;
    show_waitbar = false;

    % Define norm
    normfunc = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

    lambda = lambdas.Value(k);   

    E1_local  = cell(2,1);
    E2_local  = cell(2,1);

    local_segments = cell(n3,1);
    
    for scenario = 1:n2
        
        % Solve linear system for each segment
        for j = 1:n3
            idx = (k-1) * (n2 * n3) + (scenario-1) * n3 + j;
            local_segments{j} = compute_segments{idx}
        end
       
        %fprintf("\nIt took %.4f seconds to solve all the linear systems.\n\n",stop)

        % Compute the fields
        [~, ~,  Einc, ~, Eref, ~, Escat, ~] = compute_fields(coord_obs.Value, local_segments, far_field_approximation, scenario, lambda, show_waitbar);
        E1_local{scenario} = Einc
        E2_local{scenario} = Eref + Escat
        %E2_local{scenario} = Escat
        
        
    end

    for j = 1:length(betas.Value)
        beta = betas.Value(j);
        
    
        E1 = E1_local{1}*cos(beta) + E1_local{2}*sin(beta);
        E2 = E2_local{1}*cos(beta) + E2_local{2}*sin(beta);
        %RE_local(j,:) = E2;
        RE_local(j,:) = normfunc(abs(E1)).^2./normfunc(abs(E2));
       
    end
    RE(k,:,:) = RE_local;


end
%tocBytes(gcp)
%stop = toc;
%fprintf('It took %.4f seconds to comupute the reflectance from computed segments.\n',stop)















return
ticBytes(gcp);

betas = parallel.pool.Constant(betas);
coord_obs = parallel.pool.Constant(coord_obs);
parfor k = 1:length(lambdas.Value)
    M = length(coord_obs.Value.x);
    RE_local = zeros(length(betas.Value),M);
    RH_local = zeros(length(betas.Value),M);

    far_field_approximation = false;
    show_waitbar = false;

    % Define norm
    normfunc = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

    warning('off', 'MATLAB:singularMatrix');  % Suppresses singular matrix warnings
    warning('off','MATLAB:nearlySingularMatrix')

    lambda = lambdas.Value(k);   

    E1_local  = cell(2,1);
    E2_local  = cell(2,1);

    % Create a local copy of segments to avoid classification issues
    local_segments = segments; 

    
    for scenario = 1:2

        % Solve linear system for each segment
        for j = 1:length(segments)
            local_segments{j} = forward(local_segments{j},scenario,lambda,false);
        end
       
        %fprintf("\nIt took %.4f seconds to solve all the linear systems.\n\n",stop)

         % Compute the fields
        [~, ~,  Einc, ~, Eref, ~, Escat, ~] = compute_fields(coord_obs.Value, local_segments, far_field_approximation, scenario, lambda, show_waitbar);
        E1_local{scenario} = Einc;
        E2_local{scenario} = Eref + Escat;
  

    end

    for j = 1:length(betas.Value)
        beta = betas.Value(j);

        E1 = E1_local{1}*cos(beta) + E1_local{2}*sin(beta);
        E2 = E2_local{1}*cos(beta) + E2_local{2}*sin(beta);

        RE_local(j,:) = normfunc(abs(E1)).^2./normfunc(abs(E2));
    end
    RE(k,:,:) = RE_local;

end
tocBytes(gcp)


RE = zeros(length(lambdas),length(betas),length(coord_obs.x));
RH = zeros(length(lambdas),length(betas),length(coord_obs.x));
%wb = waitbar(0,'Computing reflectance...');
%br = 0;
return
ticBytes(gcp);
lambdas = parallel.pool.Constant(lambdas);
betas = parallel.pool.Constant(betas);
coord_obs = parallel.pool.Constant(coord_obs);
parfor k = 1:length(lambdas.Value)
    M = length(coord_obs.Value.x);
    RE_local = zeros(length(betas.Value),M);
    RH_local = zeros(length(betas.Value),M);

    far_field_approximation = false;
    show_waitbar = false;

    % Define norm
    normfunc = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

    warning('off', 'MATLAB:singularMatrix');  % Suppresses singular matrix warnings
    warning('off','MATLAB:nearlySingularMatrix')

    lambda = lambdas.Value(k);   

    E1_local  = cell(2,1);
    E2_local  = cell(2,1);

    % Create a local copy of segments to avoid classification issues
    local_segments = segments; 

    
    for scenario = 1:2

        % Solve linear system for each segment
        %tic;
        for j = 1:length(local_segments)
            local_segments{j} = forward(local_segments{j},scenario,lambda,false);
        end
        %stop = toc;
       
        %fprintf("\nIt took %.4f seconds to solve all the linear systems.\n\n",stop)

         % Compute the fields
        %tic;
        [~, ~,  Einc, ~, Eref, ~, Escat, ~] = compute_fields(coord_obs.Value, local_segments, far_field_approximation, scenario, lambda, show_waitbar);
        E1_local{scenario} = Einc;
        E2_local{scenario} = Eref + Escat;
        %stop = toc;
        %fprintf("\nIt took %.4f seconds to compute the fields.\n\n",stop)

        %br=br+1;
        %waitbar(br/(2*length(lambdas)),wb);

    end

    for j = 1:length(betas.Value)
        beta = betas.Value(j);

        E1 = E1_local{1}*cos(beta) + E1_local{2}*sin(beta);
        E2 = E2_local{1}*cos(beta) + E2_local{2}*sin(beta);

        RE_local(j,:) = normfunc(abs(E1)).^2./normfunc(abs(E2));
        %RH_local(j,:) = normfunc(abs(H1)).^2./normfunc(abs(H2));
    end
    RE(k,:,:) = RE_local;
    %RH(k,:,:) = RH_local;

end
tocBytes(gcp)
%close(wb);