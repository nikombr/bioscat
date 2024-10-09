function [RE] = compute_reflectance(protein_structure, total_x_grid_points, num_segments, coord_obs, betas, lambdas,far_field_approximation)
% beta: angle of displacement of polarisation
% lambda0: wavelength of incident plane wave in free space
% coord_obs: where we want to observe the reflectance

dir = '/Users/nikolinerehn/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/DTU/11. speciale/BioScat/';

if nargin < 7
    far_field_approximation = false; % If not specified, we do not use far field approximations
end

segment_filename = sprintf('%sData/segments_2D/%s_total_x_grid_points_%d_num_segments_%d.mat',dir,protein_structure,total_x_grid_points,num_segments);

if ~isfile(segment_filename) || strcmp(protein_structure,'backward')
    setup_settings(protein_structure, total_x_grid_points, num_segments);
end

% Get how the segments are planned out
load(segment_filename,'segments');


RE = zeros(length(lambdas),length(betas),length(coord_obs.x));
RH = zeros(length(lambdas),length(betas),length(coord_obs.x));
%wb = waitbar(0,'Computing reflectance...');
%br = 0;
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