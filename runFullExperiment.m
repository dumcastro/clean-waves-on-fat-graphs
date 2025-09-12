%% Run full experiment
clear all, clc, close all

%% Main arguments
width_sweep = {[5,5,5],[5,2.5,2.5]};
angle_sweep = {[0, pi - pi/24, pi + pi/12],[0, pi - pi/24, pi + pi/2 - pi/12];};
kappa_sweep = {0.1,0.15,0.2,0.25,0.3};

parameter_station % Go through preferred secondary arguments

for kk = 1:length(width_sweep)
    widths = width_sweep{kk};
    for ii = 1:length(kappa_sweep)
        kappa = kappa_sweep{ii};
        lambda_f = widths(1)/kappa;
        travel_distance = 25;
        Lx = lambda_f * (travel_distance + 1) / 2;
    
        for jj = 1:length(angle_sweep)
            angles = angle_sweep{jj};
    
            createFatGraph(Lx, widths, angles,graph_options);
    
            evolveWave(kappa, widths, angles,wave_options);
    
            clc, close all
        end
    
    end
end


