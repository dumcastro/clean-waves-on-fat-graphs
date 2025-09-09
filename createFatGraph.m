function [alpha] = createFatGraph(Lx, widths, angles, options)
%{
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
% Outputs:
%   fg      - FatGraph object
%   f_tilde - SC mapping object
%   C       - Canonical domain polygon
%   J       - Jacobian determinant matrix
%   w       - Grid points in canonical domain
%   z       - Grid points in physical domain
%
% Optional parameters (options struct):
%   .dxi        - Grid spacing in xi direction (default = 0.03)
%   .dzeta      - Grid spacing in zeta direction (default = 0.065)
%   .ep         - Domain extension parameter (default = 0.01)
%   .want_save  - Flag to save results (default = false)
%   .plot_flag  - Flag to generate plots (default = true)
%}

    % Set default parameter values
    if nargin < 4
        options = struct();
    end

    % Default graph parameters
    if nargin < 1 || isempty(Lx), Lx = 20; end
    if nargin < 2 || isempty(widths), widths = [1, 0.5, 0.5]; end
    if nargin < 3 || isempty(angles), angles = [0, 2*pi/3, 4*pi/3]; end

    % Default numerical parameters
    default_options = struct(...
        'dxi', 0.2, ...
        'dzeta', 0.2, ...
        'Nzeta', 22,...
        'ep', 0.01, ...
        'want_save', true, ...
        'plot_flag', false);

    % Merge user options with defaults
    option_names = fieldnames(default_options);
    for k = 1:length(option_names)
        if ~isfield(options, option_names{k})
            options.(option_names{k}) = default_options.(option_names{k});
        end
    end

    %% Create fat graph object
    fg = FatGraph(Lx, widths, angles)

    ver = fg.complex_vertices;

    %% Schwarz-Christoffel mapping
    if numel(widths) == 3
    sang = [0.5000, 1, 0.5000, 0.5000, 2.0000, 0.5000, 0.5000, 1.0000, 0.5000];
    else
    sang = [0.5000, 1, 0.5000, 0.5000, 1, 0.5000];  
    end
    P = polygon(ver);

    P_ep = outermollif(P,options.ep,angles,widths);

    %P_ep = polyedit(P_ep)
    

    % Create SC map to rectangle
    f_tilde = crrectmap(P_ep, sang);

    % Compute canonical domain
    C_tilde = evalinv(f_tilde, P_ep);

    % Calculate scaling factor
    %xi_lims = [min(real(vertex(C_tilde))), max(real(vertex(C_tilde)))];
    %zeta_lims = [min(imag(vertex(C_tilde)))+options.ep, max(imag(vertex(C_tilde)))-options.ep];
    zeta_lims = [min(imag(vertex(C_tilde))), max(imag(vertex(C_tilde)))];
    %Lxi_tilde = diff(xi_lims);
    Lzeta_tilde = diff(zeta_lims);
    %alpha1 = Lxi_tilde/(2*Lx);
    alpha = Lzeta_tilde/(widths(1)+2*options.ep); % This is the scaling factor

    %C_tilde = C_tilde -alpha*(options.ep*1i);
    


    C = (1/alpha)*C_tilde;
    %C = C -options.ep*1i - options.ep;

    %% Create computational grid
    xi_lims = [min(real(vertex(C))), max(real(vertex(C)))];
    %zeta_lims = [min(imag(vertex(C)))+options.ep, max(imag(vertex(C)))];
    zeta_lims = [min(imag(vertex(C)))+options.ep, max(imag(vertex(C)))-options.ep];
    
    options.dzeta = (zeta_lims(2)-zeta_lims(1))/(options.Nzeta-1); %dxi, dzeta are obtained from Nzeta choice
    options.dxi = options.dzeta;
    
    xi = xi_lims(1):options.dxi:xi_lims(2);
    zeta = zeta_lims(1):options.dzeta:zeta_lims(2);

    %% Dividing the domain in sectors
    node = ver(5); %physical node
    vert = evalinv(f_tilde,node)/alpha; % computing the critical node in canonical space
    targ_xi = real(vert);
    targ_zeta = imag(vert);

    th_xi = find(xi<targ_xi, 1, 'last');
    th_zeta = find(zeta<targ_zeta, 1, 'last');
    %w1 = alpha*w(:,1:th_xi); %branch 1 
    %w2 = alpha*w(1:th_zeta,th_xi:end); %branch 2
    %w3 = alpha*w(th_zeta+1:end,th_xi:end); %branch 3
    
    if numel(angles) == 3
        if ~(angles(2) + angles(3) == 2*pi) % In asymmetric case, adjust xi to end where shorter leg ends
        ble = real(vertex(C));
        ble = round(ble, 3);
        tmp = max(ble);
    
        ble(ble == tmp) = 0;
        tmp = max(ble);
        
        tmp = find(xi <= tmp, 1, 'last');
        tmp2 = tmp - th_xi;
    
        xi = xi(th_xi-tmp2:tmp);

        xi_lims= [xi(1),xi(end)];

        end
    end

    [Xi, Zeta] = meshgrid(xi, zeta);
    w = Xi + 1i * Zeta;
    z = eval(f_tilde, alpha*w);

    %% Compute Jacobian
    dz = evaldiff(f_tilde, alpha*w);
    J = (alpha^2)*abs(dz).^2;

    figure
    surf(real(w), imag(w), J);

    %% Preprocessing bugged Jacobian values

    %{
    sz = size(J);
    gap = 50;
    Nzeta = sz(1); Nxi = sz(2);

    %med1 = mean(mean(J(:,1:th_xi - gap)));
    %J(:,1:th_xi - gap) = med1*ones(Nzeta,th_xi-gap);


    


    med2 = median(median(J(1:th_zeta,th_xi+gap:tmp)));
    J(1:th_zeta,th_xi+gap:end) = med2*ones(th_zeta, Nxi-th_xi-gap+1);

    med3 = median(median(J(th_zeta:end,th_xi + gap:end)));
    J(th_zeta+1:end,th_xi + gap:end) = med3*ones(Nzeta-th_zeta, Nxi-th_xi-gap+1);
    

    %med2 = median(median(J(1:th_zeta,th_xi+gap:end)));
    %J(J == 0) = med2;

    
    figure
    surf(real(w), imag(w), J);
    %}

    %% Plot results if requested
    if options.plot_flag
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

        figure;
        surf(real(w), imag(w), J);
        xlabel('Real Axis');
        ylabel('Imaginary Axis');
        title('Jacobian Determinant of SC Transformation');
    end

    %% Save data if requested
    if options.want_save
        ang_display = round(angles, 3);
        save(['GraphData/widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat'])
    end
end
