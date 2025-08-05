%% Run full experiment
clear all, clc, close all

%% Main arguments
widths = [5, 2.5, 2.5];
angles = [0, pi - pi/24, pi + pi/24]; %angles calculated
kappa = 0.15;

lambda_f = widths(1)/kappa;
travel_distance = 3;
Lx = lambda_f * (travel_distance + 1) / 2;

%% Secondary parameters
parameter_station % Go through preferred secondary arguments

%% Generates graph data and evolve wave
% Try to process wave data; if it fails, evolve wave; if that fails, build graph
try
    % --- Try to load and process existing wave data
    fprintf('Trying to load and process existing wave data...\n');
    processWaveData(kappa, widths, angles, wave_vis_options);

catch %ME_wave
    %fprintf('Wave data not found or processing failed: %s\n', ME_wave.message);

    try
        % --- Try to evolve wave from existing graph/domain data
        fprintf('Trying to evolve wave from graph data...\n');
        evolveWave(kappa, widths, angles, wave_options);

        % After evolving, try processing again
        processWaveData(kappa, widths, angles, wave_vis_options);

    catch %ME_evolve
        %fprintf('Graph/domain data not found or evolution failed: %s\n', ME_evolve.message);

        try
            % --- Try to recreate graph/domain
            fprintf('Recreating graph/domain...\n');
            createFatGraph(Lx, widths, angles, graph_options);

            % Optional: process graph data (e.g., compute Jacobian or visualize)
            processGraphData(Lx, widths, angles, graph_vis_options);

            % Now try evolving and processing again
            evolveWave(kappa, widths, angles, wave_options);
            processWaveData(kappa, widths, angles, wave_vis_options);

        catch %ME_create
            %fprintf('Failed to create graph/domain: %s\n', ME_create.message);
        end

    end
end
