% Define circle parameters
x_c = 50;   % x-coordinate of the circle center
y_c = 50;   % y-coordinate of the circle center
r = 20;     % radius of the circle

% Define point on the circle
x_p = x_c + r * cos(pi/4);   % x-coordinate of the point (example point at 45 degrees)
y_p = y_c + r * sin(pi/4);   % y-coordinate of the point

% Define the boundary function y = f(x)
f = @(x) 0.01 * (x - 50).^2 + 80;   % Example: parabola shifted up and centered at x = 50

% Compute the outward normal vector (radial direction)
normal_x = (x_p - x_c) / r;
normal_y = (y_p - y_c) / r;

% Define the parametric line equation: x = x_p + t * normal_x, y = y_p + t * normal_y
% We want to find t such that y_p + t * normal_y = f(x_p + t * normal_x)

% Define the function to solve for t
g = @(t) y_p + t * normal_y - f(x_p + t * normal_x);

% Use fzero to find the root of g(t)
t_intersection = fzero(g, 1);  % Initial guess for t

% Compute the intersection point on the boundary
x_b = x_p + t_intersection * normal_x;
y_b = y_p + t_intersection * normal_y;

% Display the circle, boundary, and projection
figure;
hold on;

% Plot the boundary function
fplot(f, [0, 100], 'k-', 'LineWidth', 1.5);  % Plot y = f(x)
hold on
% Draw the circle
viscircles([x_c, y_c], r, 'LineStyle', '--');

% Mark the point on the circle
plot(x_p, y_p, 'ro');

% Plot the line from the circle to the boundary
plot([x_p, x_b], [y_p, y_b], 'b-');

% Mark the intersection point
plot(x_b, y_b, 'go');

axis equal;
title('Projection of Circle Point to Curved Boundary');
hold off;
