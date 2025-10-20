%% Run full experiment (parameter sweep)
clear all, clc, close all

addpath('External/')

%% Sweep parameters
widthTriplets = {[5,5,6],[5,5,7],[5,5,8],[5,5,9]};

theta2 = pi/6;

kappas = 0.1:0.1:0.4;
thetas3 = pi/6;

parameterStation

%% Construct data

parameterSweep(widthTriplets, kappas, theta2, thetas3)

%% Visualize data

%parameterSweepVis(widthTriplets, kappas, theta2, thetas3, ...
    %travelDistance, colorGrid, deltaHeightPlot)








