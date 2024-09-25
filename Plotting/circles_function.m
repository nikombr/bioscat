

clear; clc; close all;

% Parameters for the circles (1D projections of their centers and radii)
circles = [1, 0.5;  % [x1, r1]
           3, 1.0;  % [x2, r2]
           6, 0.75;  % [x3, r3]
           8,    1]; % [x4, r4]
circles = [3, 1.0;
           8,    1]; % [x4, r4]
nCircles = size(circles, 1);

% Define the x-axis over which the function will be computed
x_range = linspace(0,10,200);

% Initialize the smooth function
smoothFunction = zeros(size(x_range));


% For each circle, add a Gaussian function representing the 1D interval
for i = 1:nCircles
    x_center = circles(i, 1);
    radius = circles(i, 2);

    % Compute the distance of each x point to the circle's center
        distanceFromCenter = abs(x_range - x_center);

        % Determine which x points are within the circle's width (|x - x_i| <= r_i)
        insideCircle = distanceFromCenter <= radius;

        % Compute the corresponding y-values for points inside the circle's projection
        % Use the equation y = sqrt(r^2 - (x - x_center)^2)
        y_values = radius + sqrt(radius^2 - (x_range - x_center).^2);

        % Set y = 0 for points outside the circle's projection
        y_values(~insideCircle) = 0;

        % Accumulate the smooth function for all circles
        smoothFunction = smoothFunction + y_values;
    
    % Gaussian function centered at x_center with width proportional to radius
    %smoothFunction = smoothFunction + exp(-((x_range - x_center).^2) / (2 * (radius^2)));
end

% Normalize the function (optional)
%smoothFunction = smoothFunction / max(smoothFunction);

for i = 1:10
smoothFunction = smooth(smoothFunction,5);
end


% Plot the result
figure;
plot(x_range, smoothFunction, 'LineWidth', 2);
title('Smooth Function Over Projected Circles in 1D');
xlabel('x');
ylabel('f(x)');
grid on;
axis equal
hold on


for i = 1:nCircles
    circle = circles(i,:);
    x = circle(1);
    r = circle(2);
    y = r;
    fplot(@(t) r*sin(t)+x, @(t) r*cos(t)+y);
end

xint = linspace(0,10,150);
yint = interp1(x_range,smoothFunction,xint);

for i = 1:length(xint)
    if i == length(xint)
        fprintf("(%.4f,%.4f);",xint(i),yint(i))
    else
        fprintf("(%.4f,%.4f) -- ",xint(i),yint(i))
    end
    if mod(i,5) == 0
        fprintf("\n")
    end
end


