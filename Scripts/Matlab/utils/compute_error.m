function error = compute_error(true_field,predicted_field,relative)


% Reshape fields
true_field      = reshape(true_field,     [],3);
predicted_field = reshape(predicted_field,[],3);

% Compute difference
diff = abs(true_field - predicted_field);

% Define norm
normfunc = @(x) sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);

% Compute absolute error
error = normfunc(diff);

% Normalize if we want to compute the relative error
if relative
    error = error./normfunc(abs(true_field));
end





