function [H,U,V,h] = evolveWave2(kappa, Lx, widths, angles, options)
%EVOLVEWAVE Linear wave evolution
%   Choose geometry of fat graph and solve wave evolution in that domain

% Inputs:
%   kappa   - width-to-wavelenght regime 
%   Lx      - Length of each branch
%   widths  - Array of branch widths [main, branch1, branch2] 
%   angles  - Array of branch angles in radians (default = [0, 2*pi/3, 4*pi/3])
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
    if nargin < 1 || isempty(kappa), kappa = 0.05; end
    if nargin < 2 || isempty(widths), widths = [5, 2.5, 2.5]; end
    if nargin < 3 || isempty(angles), angles = [0, 2*pi/3, 4*pi/3]; end

    % Default numerical parameters
    default_options = struct(...
        'dt', 0.02,...
        'T', 200,...
        'travel_distance', 20,...
        'frames', 25,...
        'want_save', true,...
        'point_source', false);

    % Merge user options with defaults
    option_names = fieldnames(default_options);
    for k = 1:length(option_names)
        if ~isfield(options, option_names{k})
            options.(option_names{k}) = default_options.(option_names{k});
        end
    end

%% Load pre-generated graph data
[~, graphName] = standardNaming(Lx, widths, angles, kappa);

load(graphName,'J','xi', 'zeta', 'Xi', 'Zeta', 'dxi','dzeta');
    
%% Setting parameters and initial data
dt = options.dt;
T = options.T;

% Define Gaussian pulse parameters
lambda_f = widths(1)/kappa;
sigma = lambda_f / 6.065; % (this is what sigma must be to make the effective wavelength in canonical
% space to be alpha*lambda_f)

xi0 = ((xi(end))+(xi(1)))/ 2;
zeta0 = ((zeta(end))+(zeta(1)))/ 2;

a = 0.1; %pulse height

h = a*exp(-(Xi-xi0/2).^2/ (2 * sigma^2));
u = zeros(size(h));        % Initial velocity
v = h.*J.^(1/2); % necessary velocity for unidirectional solution (right-going mode only)

if options.point_source
    sigma = 0.3;
    xic = xi0/2; zetac = zeta0/2;
    h  = a*exp((-((Xi-xic).^2) - (Zeta-zetac).^2)/sigma^2);
    v = u;
end

h = neumann_correction(h);
u = impermiability_u(u);
v = impermiability_v(v);

%% Setting wave data

H = reshape(h, numel(J),1);
U = reshape(u, numel(J),1);
V = reshape(v, numel(J),1);

%% Main loop % RK4 time-stepping
%dist = 0;
t = 0; iter=0;
t_array = 0:dt:T;
tol = 0.001*a;
%while t < T
while max(h(:, end)) < tol
    k1_u = -dt * (circshift(h, [ -1 0]) - circshift(h, [ 1 0])) / (2 * dxi);
    k1_v = -dt * (circshift(h, [ 0 -1]) - circshift(h, [0 1])) / (2 * dzeta);
    k1_h = -dt * ((circshift(u, [ -1 0]) - circshift(u, [ 1 0])) / (2 * dxi) +...
        (circshift(v, [ 0 -1]) - circshift(v, [ 0 1])) / (2 * dzeta))./J;
    
    h1 = h + 0.5 * k1_h;
    u1 = u + 0.5 * k1_u;
    v1 = v + 0.5 * k1_v;
    
    h1 = neumann_correction(h1);
    u1 = impermiability_u(u1);
    v1 = impermiability_v(v1);
    
    k2_u = -dt * (circshift(h1, [ -1 0]) - circshift(h1, [ 1 0])) / (2 * dxi);
    k2_v = -dt * (circshift(h1, [ 0 -1]) - circshift(h1, [ 0 1])) / (2 * dzeta);
    k2_h = -dt * ((circshift(u1, [ -1 0]) - circshift(u1, [ 1 0])) / (2 * dxi) + ...
        (circshift(v1, [ 0 -1]) - circshift(v1, [ 0 1])) / (2 * dzeta))./J;
    
    h2 = h + 0.5 * k2_h;
    u2 = u + 0.5 * k2_u;
    v2 = v + 0.5 * k2_v;
    
    h2 = neumann_correction(h2);
    u2 = impermiability_u(u2);
    v2 = impermiability_v(v2);
    
    k3_u = -dt * (circshift(h2, [ -1 0]) - circshift(h2, [ 1 0])) / (2 * dxi);
    k3_v = -dt * (circshift(h2, [ 0 -1]) - circshift(h2, [ 0 1])) / (2 * dzeta);
    k3_h = -dt * ((circshift(u2, [ -1 0]) - circshift(u2, [ 1 0])) / (2 * dxi) + ...
        (circshift(v2, [ 0 -1]) - circshift(v2, [ 0 1])) / (2 * dzeta))./J;
   
    h3 = h + k3_h;
    u3 = u + k3_u;
    v3 = v + k3_v;

    h3 = neumann_correction(h3);
    u3 = impermiability_u(u3);
    v3 = impermiability_v(v3);
    
    k4_u = -dt * (circshift(h3, [ -1 0]) - circshift(h3, [ 1 0])) / (2 * dxi);
    k4_v = -dt * (circshift(h3, [0 -1]) - circshift(h3, [ 0 1])) / (2 * dzeta);
    k4_h = -dt * ((circshift(u3, [ -1 0]) - circshift(u3, [ 1 0])) / (2 * dxi) + ...
        (circshift(v3, [ 0 -1]) - circshift(v3, [ 0 1])) / (2 * dzeta))./J;
    
    u = u + (1/6) * (k1_u + 2 * k2_u + 2 * k3_u + k4_u);
    v = v + (1/6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
    h = h + (1/6) * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);
    
    v = impermiability_v(v);
    h = neumann_correction(h);
    u = impermiability_u(u);
    
    %% Saves selected number of frames
    if mod(iter,floor(numel(t_array)/options.frames))==0 || t == T-1-dt
        
        H = [H, reshape(h,numel(J),1)];
        %U = [U, reshape(u,numel(J),1)];
        %V = [V, reshape(v,numel(J),1)];
                
    end
    
    % Update time
    t = t + dt;
    hh = h(floor(end/2),:);   
    [~,loc] = findpeaks(hh,'MinPeakHeight', 0.7*a); % finds index for which final wave prof peaks (center of final gaussian)
    x0f = xi(loc);
    dist = abs(xi0 - x0f); %gets distance between current and initial peak positions
    
    iter=iter+1;
    max(h(:, end))
    
    if t > T
        message = ['Exceeded final time ceiling. T = ', num2str(T)]
        break
    end
end

%% Save data if requested
    if options.want_save

       waveName = standardNaming(Lx, widths, angles, kappa);
       save(waveName,'H','U', 'V','h')

    end
  
% Aux functions    
%% Exterior boundary contidtion functions
function v = impermiability_v(v)
    % Set v to zero along the boundaries
    v(:, 1) = 0; % Left boundary
    v(:, end) = 0; % Right boundary
    
    %v(th_zeta+1,th_xi-1:end)=0;
    %v(th_zeta,th_xi-1:end)=0;
    
end

function u = impermiability_u(u)
    % Set v to zero along the boundaries
    u(1, :) = 0; % Bottom boundary
    u(end, :) = 0; % Top boundary
end

function h=neumann_correction(h)

    h(:, 1) = h(:, 2); % Left boundary
    h(:, end) = h(:, end-1); % Right boundary   
    
    h(1, :) = h(2, :); % Bottom boundary
    h(end, :) = h(end-1, :); % Top boundary
end

end

