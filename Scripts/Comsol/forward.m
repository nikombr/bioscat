close all force
clear all
clc
beep off
format compact
format long

% Load nanostructure
xs=load("../../Data/nanostructures/demoleus2x2_2D_x_300.txt");
fs=load("../../Data/nanostructures/demoleus2x2_2D_f_300.txt");

% Load modules
import com.comsol.model.*
import com.comsol.model.util.*

figure, plot(xs,fs), grid on, axis equal
surface=[xs' fs'];
save('surface.txt','surface','-ascii')

% Load Comsol model
model = mphload('forward.mph');

% Solve
model.sol('sol1').runAll;

% Get solution
Ex = mpheval(model,'emw.Ez')