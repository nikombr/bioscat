
clear; close all; clc;

filename = "demoleus2x2"; % Retinin2x2 demoleus2x2

load(sprintf("../../Data/%s_2D.mat",filename))


m = length(X);
num_segments = 10;
n = round(m/num_segments);
n = 50;
m = num_segments * n;

x = linspace(min(X),max(X),m+1);
d = x(2)-x(1);

r = (max(x) - min(x))/(2*num_segments);
alpha = 0.6;

y = interp1(X,Y,x);

figure('Renderer', 'painters', 'Position', [400 400 800 250]);
plot(x,y,'LineWidth',1.5)
xlim([min(x),max(x)])
grid on
hold on

segments = cell(num_segments,1);

figure('Renderer', 'painters', 'Position', [400 400 1000 800]);
for k = 1:num_segments
    
    segx = x((k-1)*n+1:k*n+1);
    segy = y((k-1)*n+1:k*n+1);
    
    segments{k} = compute_auxiliary_sources_and_test_points_2D(segx,segy,alpha);

    plot(segments{k}.x,segments{k}.y,'r.','LineWidth',1.5,'Markersize',10)
    hold on
    
    plot(segments{k}.x_int,segments{k}.y_int,'b.','LineWidth',1.5,'Markersize',10)
    
    plot(segments{k}.x_ext,segments{k}.y_ext,'go','LineWidth',1.5,'Markersize',5)
end
%xlim([min(x_segment),max(x_segment)])
xlim([-0.5*10^(-7),20.5*10^(-7)])
%ylim([-0.5*10^(-7),5*10^(-7)])
ylim([-0.5*10^(-7),5*10^(-7)])
grid on
axis equal


