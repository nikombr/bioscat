close all force
clear all
clc
beep off
format compact
format long

xs=load("demoleus2x2_2D_x_300.txt");
fs=load("demoleus2x2_2D_f_300.txt");

figure, plot(xs,fs), grid on, axis equal
surface=[xs' fs'];
save('surface.txt','surface','-ascii')