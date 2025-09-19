%% THIS IS THE MAIN/DEMO SCRIPT
clear all, clc, close all

addpath('External/')

% Main arguments
widths = [5, 5, 5];

theta2 = pi/30; theta3 = 5*pi/12;
angles = [0, pi - theta2, pi + theta3];

kappa = 0.25;

% Secondary parameters
parameter_station % Go through preferred secondary arguments
lambda_f = widths(1)/kappa;
Lx = lambda_f * (travel_distance + 1) / 2;
%Lx = 200;

%% Testing create fat graph
%createFatGraph(Lx, widths, angles,graph_options);

%% Testing process Graph data
%processGraphData(Lx, widths, angles,graph_vis_options)

%% Testing evolveWave
%
if numel(widths) == 3
evolveWave(kappa, Lx, widths, angles,wave_options)
elseif numel(widths) == 2
evolveWave2(kappa, Lx, widths, angles,wave_options)   
end
%}
%% Testing processWave
processWaveData(kappa, Lx, widths, angles,wave_vis_options)