% Copyright (C) 2025 Eduardo Castro

% Author: Eduardo Castro <eduardocastro@MacBook-Air-de-Eduardo.local>
% Created: 2025-06-03

function [] = processGraphData(Lx, widths, angles, options)
%CREATEFATGRAPH Creates a fat graph and computes Schwarz-Christoffel mapping
%   [fg, f_tilde, C, J, w, z] = createFatGraph(Lx, widths, angles) creates a Y-shaped
%   fat graph with specified dimensions and computes its SC transformation.
%
%   [fg, f_tilde, C, J, w, z] = createFatGraph(Lx, widths, angles, options) allows
%   customization of numerical parameters.
%
% Inputs:
%   Lx      - Length of main branch (default = 10)
%   widths  - Array of branch widths [main, branch1, branch2] (default = [1, 0.4, 0.6])
%   angles  - Array of branch angles in radians (default = [0, 2*pi/3, 4*pi/3])
%   options - Structure with optional parameters (see below)
%
% Optional parameters (options struct):
%   .ep         - Domain extension parameter (default = 0.01)
%   .want_save  - Flag to save results (default = false)
%   .plot_flag  - Flag to generate plots (default = true)

    % Set default parameter values
    if nargin < 4
        options = struct();
    end

    % Default graph parameters
    if nargin < 1 || isempty(Lx), Lx = 10; end
    if nargin < 2 || isempty(widths), widths = [1, 0.5, 0.5]; end
    if nargin < 3 || isempty(angles), angles = [0, 2*pi/3, 4*pi/3]; end

    % Default numerical parameters
    default_options = struct(...
        'ep', 0.01,...
        'jmp_xi', 1,...
        'jmp_zeta', 1,...
        'plotJ', true,...
        'plotMap', false,...
        'visualizeGrid', true);

    % Merge user options with defaults
    option_names = fieldnames(default_options);
    for k = 1:length(option_names)
        if ~isfield(options, option_names{k})
            options.(option_names{k}) = default_options.(option_names{k});
        end
    end

    %% Load fat Graph
    ang_display = round(angles, 3);
    %data = load(['GraphData/widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat']);
    load(['GraphData/widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat'],...
                'w', 'J', 'z','th_xi','th_zeta','xi', 'zeta', 'Xi', 'Zeta', 'dxi','dzeta', 'vert');

    %% Plot results if requested
    if options.plotMap
        figure;
        subplot(1, 3, 1);
        plot(P, 'b', 'LineWidth', 2, 'k');
        title('Physical Region');

        subplot(1, 3, 2);
        plot(C_tilde, 'r', 'LineWidth', 2);
        title('Numerical canonical domain');

        subplot(1, 3, 3);
        plot(C, 'y', 'LineWidth', 2);
        title('Canonical domain width = 1');
    end
    
    
    if options.plotJ
        figure;
        
        tmp_Xi = real(w);
        tmp_Zeta = imag( w);
        J = J;
        
        jmp_xi = options.jmp_xi;
        jmp_zeta = options.jmp_zeta;
        
        tmp_Xi = tmp_Xi(1:jmp_zeta:end,1:jmp_xi:end);
        tmp_Zeta = tmp_Zeta(1:jmp_zeta:end,1:jmp_xi:end);
        
        %surf(real(w), imag(w), J);
        surf(tmp_Xi,tmp_Zeta,J(1:jmp_zeta:end,1:jmp_xi:end));
        xlabel('\xi'); ylabel('\zeta', 'Rotation', 0);
        xlim([real(vert)-10,real(vert)+10])
        %title('Jacobian Determinant of SC Transformation');

        set(gca, 'FontSize', 16)   % makes axis numbers larger
    end
    

%% Visualize Grid
    if options.visualizeGrid
        % Visualize with scatter
        gzoom = .5;
        
        figure;
        scatter(Xi(:), Zeta(:), 100, 'filled');  % Flatten grids and plot
        xlim([real(vert)-gzoom,real(vert)+gzoom])
        ylim([imag(vert)-gzoom,imag(vert)+gzoom])
        %title('Meshgrid Visualization');
        xlabel('\xi'); ylabel('\zeta', 'Rotation', 0);
        grid on; hold on,
        
        xcoord = [real(vert), xi(end)];
        ycoord = [imag(vert), imag(vert)];

        scatter(real(vert), imag(vert) ,120,'red','d','filled')

        plot(xcoord,ycoord, 'r-', 'LineWidth', 2)

        set(gca, 'FontSize', 16)   % makes axis numbers larger

        


    end
end

