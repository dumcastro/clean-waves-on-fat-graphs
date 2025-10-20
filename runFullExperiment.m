%% Run full experiment (parameter sweep)
clear all, clc, close all

addpath('External/')

%% Sweep parameters
widthTriplets = {[5,5,3],[5,5,5],[5,5,7]};
%widthTriplets = {[5,5,5]};

theta2 = pi/30;

kappas = 0.1:0.2:0.5;
%kappas = 0.3;
%thetas3 = pi/12:pi/6:5*pi/12;

thetas3 = 5*pi/12;


parameterStation

%% Construct data

%parameterSweep(widthTriplets, kappas, theta2, thetas3)

%% Visualize data

parameterSweepVis(widthTriplets, kappas, theta2, thetas3, ...
    travelDistance, sweepVisOptions)








