%% Run full experiment (parameter sweep)
clear all, clc, close all

%% Main arguments
width_sweep = {[5,5,5]};

theta2 = pi/30;

anglesStart = [0, pi - theta2, pi + 0];
angles = anglesStart;

kappas = 0.25;
thetas3 = 5*pi/12;

%angle_sweep = {[0, pi - pi/30, pi + pi/10],[0, pi - pi/30, pi + pi/2 - pi/30];};
%kappa_sweep = {0.08, 0.25, 0.35};

parameter_station % Go through preferred secondary arguments

for kk = 1:length(width_sweep)
    widths = width_sweep{kk};
    for ii = 1:length(kappas)
        kappa = kappas(ii);
        lambda_f = widths(1)/kappa;
        Lx = lambda_f * (travel_distance + 1) / 2;
    
        for jj = 1:length(thetas3)
            angles(3) = anglesStart(3) + thetas3(jj);
    
            createFatGraph(Lx, widths, angles,graph_options);
    
            evolveWave(kappa, Lx, widths, angles,wave_options);
    
            close all
        end
        angles = anglesStart;

    end
end


