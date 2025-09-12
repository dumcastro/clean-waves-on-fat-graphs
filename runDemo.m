%% THIS IS THE MAIN/DEMO SCRIPT
clear all, clc, close all

addpath('External/')

% Main arguments
widths = [5, 5 ,5];
angles = [0, pi - pi/12, pi + pi/2 - pi/12];

kappa = 0.18;

lambda_f = widths(1)/kappa;
travel_distance = 4;
Lx = lambda_f * (travel_distance + 1) / 2;
%Lx = 200;

% Secondary parameters
parameter_station % Go through preferred secondary arguments

%% Testing FatGraph
%fg = FatGraph(Lx, widths, angles)

%% Testing create fat graph
%createFatGraph(Lx, widths, angles,graph_options);

%% Testing process Graph data

%processGraphData(Lx, widths, angles,graph_vis_options)

%% Testing evolveWave
%
if numel(widths) == 3
evolveWave(kappa, widths, angles,wave_options)
elseif numel(widths) == 2
evolveWave2(kappa, widths, angles,wave_options)   
end
%}
%% Testing processWave
processWaveData(kappa, widths, angles,wave_vis_options)
