%% THIS IS THE MAIN/DEMO SCRIPT
clear all, clc, close all

addpath('External/')

% Main arguments
widths = [5, 5, 5];

theta2 = pi/12; theta3 = 5*pi/12;
angles = [0, pi - theta2, pi + theta3];

kappa = 0.2;

% Secondary parameters
parameterStation % Go through preferred secondary arguments
lambda_f = widths(1)/kappa;
Lx = lambda_f * (travelDistance + 1) / 2;
%Lx = 200;

%% Testing create fat graph
%createFatGraph(Lx, widths, angles,graph_options);

%% Testing process Graph data
processGraphData(Lx, widths, angles,graph_vis_options)

%% Testing evolveWave
%evolveWave(kappa, Lx, widths, angles,wave_options);

%% Testing processWave
%processWaveData(kappa, Lx, widths, angles,wave_vis_options)