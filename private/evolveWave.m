function [H,U,V,h] = evolveWave(kappa, Lx, widths, angles, options)
%EVOLVEWAVE Linear wave evolution
%   Choose geometry of fat graph and solve wave evolution in that domain

% Inputs:
%   kappa   - width-to-wavelenght regime
%   Lx      - Length of each branch, all equal for now
%   widths  - Array of branch widths [main, branch2, branch3] 
%   angles  - Array of branch angles in radians
%   options - Structure with optional parameters (see below)

%
% Optional parameters (options struct):
%   .dt         - Time step size (default = 0.02) 
%   .ep         - Domain extension parameter (default = 0.01)
%   .want_save  - Flag to save results (default = true)
%   .T          - Final time of experiment (default = 115)
%   .frames     - How many frames will be saved (default = 25)

    % Set default parameter values
    if nargin < 4
        options = struct();
    end

    % Default graph parameters
    if nargin < 1 || isempty(kappa), kappa = 0.1; end
    if nargin < 2 || isempty(widths), widths = [5, 5, 5]; end
    if nargin < 3 || isempty(angles), angles = [0, 2*pi/3, 4*pi/3]; end

    % Default numerical parameters
    default_options = struct(...
        'dt', 0.02,...
        'T', 200,...
        'travel_distance', 20,...
        'frames', 25,...
        'want_save', true,...
        'point_source', false,...
        'transverse_wave', false);

    % Merge user options with defaults
    option_names = fieldnames(default_options);
    for k = 1:length(option_names)
        if ~isfield(options, option_names{k})
            options.(option_names{k}) = default_options.(option_names{k});
        end
    end

%% Load pre-generated graph data
ang_display = round(angles, 3);
load(['GraphData/widths=', mat2str(widths), '_angles=', mat2str(ang_display),'_length=',...
            mat2str(Lx), '.mat'],...
            'w', 'J', 'z','th_xi','th_zeta','xi', 'zeta', 'Xi', 'Zeta', 'dxi','dzeta');
    
%% Setting parameters and initial data
dt = options.dt;
T = options.T;

% Define Gaussian pulse parameters
lambda_f = widths(1)/kappa;
comp_efetivo_can = lambda_f;
sigma = comp_efetivo_can / 6.065; % (this is what sigma must be to make the effective wavelength in canonical
% space to be alpha*lambda_f)

xi0 = ((xi(end))+(xi(1)))/ 2;
zeta0 = ((zeta(end))+(zeta(1)))/ 2;

a = 0.1; %pulse height

h = a*exp(-(Xi-xi(floor(numel(xi)/4))).^2/ (2 * sigma^2));
u = zeros(size(h));        % Initial velocity
v = h.*J.^(1/2); % necessary velocity for unidirectional solution (right-going mode only)

if options.point_source
    sigma = 0.4;
    xic = 1.1*xi0; zetac = 1.4*zeta0;
    h  = a*exp((-((Xi-xic).^2) - (Zeta-zetac).^2)/sigma^2);
    v = u;
end

if options.transverse_wave %(for debug purposes)
    % Initial data:
    h = a*exp(-(Zeta-4).^2/ (2* kappa^4 * sigma^2));
    u = -h; % necessary velocity for unidirectional solution (down-going mode only)
    v = zeros(size(h));
end

[h, u, v] = enforce_barrier(h, u, v, th_zeta, th_xi);  % Enforce slit barrier

%% Setting wave data
H = reshape(h, numel(J),1);
U = reshape(u, numel(J),1);
V = reshape(v, numel(J),1);

%% Main loop % RK4 time-stepping
t = 0; iter = 0;
t_array = 0:dt:T;
tol = 0.0001*a;
%while t < T
while max(h(:, end)) < tol
    k1_u = -dt * (circshift(h, [ -1 0]) - circshift(h, [ 1 0])) / (2 * dxi);
    k1_v = -dt * (circshift(h, [ 0 -1]) - circshift(h, [0 1])) / (2 * dzeta);
    k1_h = -dt * ((circshift(u, [ -1 0]) - circshift(u, [ 1 0])) / (2 * dxi) +...
        (circshift(v, [ 0 -1]) - circshift(v, [ 0 1])) / (2 * dzeta))./J;
    
    h1 = h + 0.5 * k1_h;
    u1 = u + 0.5 * k1_u;
    v1 = v + 0.5 * k1_v;
    
    [h1, u1, v1] = enforce_barrier(h1, u1, v1, th_zeta, th_xi);  % Enforce slit barrier on partial rk4 sums
    
    k2_u = -dt * (circshift(h1, [ -1 0]) - circshift(h1, [ 1 0])) / (2 * dxi);
    k2_v = -dt * (circshift(h1, [ 0 -1]) - circshift(h1, [ 0 1])) / (2 * dzeta);
    k2_h = -dt * ((circshift(u1, [ -1 0]) - circshift(u1, [ 1 0])) / (2 * dxi) + ...
        (circshift(v1, [ 0 -1]) - circshift(v1, [ 0 1])) / (2 * dzeta))./J;
    
    h2 = h + 0.5 * k2_h;
    u2 = u + 0.5 * k2_u;
    v2 = v + 0.5 * k2_v;
    
    [h2, u2, v2] = enforce_barrier(h2, u2, v2, th_zeta, th_xi);  % Enforce slit barrier on partial rk4 sums
    
    k3_u = -dt * (circshift(h2, [ -1 0]) - circshift(h2, [ 1 0])) / (2 * dxi);
    k3_v = -dt * (circshift(h2, [ 0 -1]) - circshift(h2, [ 0 1])) / (2 * dzeta);
    k3_h = -dt * ((circshift(u2, [ -1 0]) - circshift(u2, [ 1 0])) / (2 * dxi) + ...
        (circshift(v2, [ 0 -1]) - circshift(v2, [ 0 1])) / (2 * dzeta))./J;
    
    h3 = h + k3_h;
    u3 = u + k3_u;
    v3 = v + k3_v;
    
    [h3, u3, v3] = enforce_barrier(h3, u3, v3, th_zeta, th_xi);  % Enforce slit barrier on partial rk4 sums
    
    k4_u = -dt * (circshift(h3, [ -1 0]) - circshift(h3, [ 1 0])) / (2 * dxi);
    k4_v = -dt * (circshift(h3, [0 -1]) - circshift(h3, [ 0 1])) / (2 * dzeta);
    k4_h = -dt * ((circshift(u3, [ -1 0]) - circshift(u3, [ 1 0])) / (2 * dxi) + ...
        (circshift(v3, [ 0 -1]) - circshift(v3, [ 0 1])) / (2 * dzeta))./J;
    
    u = u + (1/6) * (k1_u + 2 * k2_u + 2 * k3_u + k4_u);
    v = v + (1/6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    h = h + (1/6) * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);
    
    [h, u, v] = enforce_barrier(h, u, v, th_zeta, th_xi);  % Enforce slit barrier
    
    %% Saves selected number of frames
    if mod(iter,floor(numel(t_array)/options.frames))==0 || t == T-1-dt
        
        H = [H, reshape(h,numel(J),1)];
        U = [U, reshape(u,numel(J),1)];
        V = [V, reshape(v,numel(J),1)];
        
    end
    
    % Update time
    t = t + dt;
    
    iter=iter+1;
    
    disp(['W. height at the end: ', num2str(max(h(:, end))), ' out of', num2str(tol)])
    disp(['Time: ', num2str(t)])
    
    if t > T
        disp(['Exceeded final time ceiling. T = ', num2str(T)])
        break
    end

    if max(h(:, end)) >= tol
        disp('Wave front reaching the end of comp. domain. Stop.')
        break
    end

end

%% Save data if requested
    if options.want_save
        ang_display = round(angles, 3);
        save(['WaveData/kappa=', num2str(kappa),'_widths=', mat2str(widths), '_angles=', mat2str(ang_display), '_length=',...
            mat2str(Lx), '.mat'], ...
            'H','U', 'V','h');
    end
       
%% Boundary condition function
function [h,u, v] = enforce_barrier(h, u, v, b1, b2)
    % Slit barrier enforcement
    h(b1+1,b2:end) = h(b1+2,b2:end);   
    u(b1+1,b2:end) = 0;
    %v(b1+1,b2:end) = 0;
    
    h(b1-1, b2:end) = h(b1, b2:end); 
    u(b1,b2:end) = 0;
    %v(b1,b2:end) = 0;

    % Outer walls conditions
    v(:, 1) = 0; % Left boundary
    v(:, end) = 0; % Right boundary

    u(1, :) = 0; % Bottom boundary
    u(end, :) = 0; % Top boundary

    h(:, 1) = h(:, 2); % Left boundary
    h(:, end) = h(:, end-1); % Right boundary   
    
    h(1, :) = h(2, :); % Bottom boundary
    h(end, :) = h(end-1, :); % Top boundary
    %}
    
end

end

