
% Run "comsolmatlab" in the terminal to start matlab session with access to comsol

clear; close all; clc;

n = 4
t = n + 2
"hej"

% Set number of threads
ont = maxNumCompThreads('automatic');
nnt = maxNumCompThreads();
fprintf('Switching from %d to %d threads!\n', ont, nnt);

% Load workflow params
load params