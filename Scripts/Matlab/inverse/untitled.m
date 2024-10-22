clear; close all; clc;

u = pi;
v = pi/2+pi/3;
x = [cos(u)*cos(v) sin(u)*cos(v) sin(v)];

figure;
quiver3(0,0,0,x(1),x(2),x(3))
hold on
quiver3(0,0,0,-1,0,0)
xlabel('x')

ylabel('y')
zlabel('z')