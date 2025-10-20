%% Run full experiment (parameter sweep)
clear all, clc, close all

addpath('External/') %make sure dependencies are in MATLAB's path

%% Sweep parameters
widthTriplets = {[5,5,3],[5,5,4],[5,5,5],[5,5,6],[5,5,7]};

kappas = 0.05:0.05:0.5;

theta2 = pi/30; % we let this leg be fixed always
thetas3 = pi/12:pi/12:7*pi/12; % this is the leg that sweeps

parameterStation % check this script

%% Construct data

%parameterSweep(widthTriplets, kappas, theta2, thetas3)

%% Visualize data

parameterSweepVis(widthTriplets, kappas, theta2, thetas3, ...
    travelDistance, sweepVisOptions)








