%% Run full experiment
clear all, clc, close all

%% Main arguments
width_sweep = {[5,5,5]};

anglesStart = [0, pi - pi/30, pi + 0];
angles = anglesStart;

kappas = 0.1:0.1:0.4;
thetas = 3*pi/12: pi/6 :pi/2-pi/12;

%angle_sweep = {[0, pi - pi/30, pi + pi/10],[0, pi - pi/30, pi + pi/2 - pi/30];};
%kappa_sweep = {0.08, 0.25, 0.35};

travel_distance = 10;

for kk = 1:length(width_sweep)
    widths = width_sweep{kk};
    parameter_station % Go through preferred secondary arguments
    for ii = 1:length(kappas)
        kappa = kappas(ii);
        lambda_f = widths(1)/kappa;
        Lx = lambda_f * (travel_distance + 1) / 2;
    
        for jj = 1:length(thetas)
            angles(3) = anglesStart(3) + thetas(jj);
    
            createFatGraph(Lx, widths, angles,graph_options);
    
            evolveWave(kappa, widths, angles,wave_options);
    
            clc, close all
        end
        angles = anglesStart;

    end
end


