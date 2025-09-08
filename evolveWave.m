function [] = evolveWave(kappa, widths, angles, options)
%EVOLVEWAVE Linear wave evolution
%   Choose geometry of fat graph and solve wave evolution in that domain

% Inputs:
%   kappa      - width-to-wavelenght regime (default = 0.05)
%   widths  - Array of branch widths [main, branch1, branch2] (default = [5, 2.5, 2.5])
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
data = load(['GraphData/widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat']);
    
%% Setting parameters and initial data
dt = options.dt;
T = options.T;
dxi = data.options.dxi;
dzeta = data.options.dzeta;
%alpha = data.alpha;
z = data.z;
w = data.w;
xi_lims = data.xi_lims;
Xi = data.Xi;
xi = data.xi;
Zeta = data.Zeta;
J = data.J;
th_zeta = data.th_zeta;
th_xi = data.th_xi-data.tmp2;

% Define Gaussian pulse parameters
lambda_f = widths(1)/kappa;
comp_efetivo_can = lambda_f;
sigma = comp_efetivo_can / 6.065; % (this is what sigma must be to make the effective wavelength in canonical
% space to be alpha*lambda_f)

zeta_lims = [imag(w(end,1)), imag(w(1,1))];


%xi0 = ((xi_lims(2))+(xi_lims(1)))/ 2;
zeta0 = ((zeta_lims(2))+(zeta_lims(1)))/ 2;

a = 0.1; %pulse height

h = a*exp(-(Xi-xi(floor(numel(xi)/4))).^2/ (2 * sigma^2));
u = zeros(size(h));        % Initial velocity
v = h.*J.^(1/2); % necessary velocity for unidirectional solution (right-going mode only)

if options.point_source
    h  = 20*a*exp(-((Xi-2*xi0+xi0/4).^2)/40 - (Zeta-zeta0/3).^2)/(20 * sigma^2);
    v = u;
end

if options.transverse_wave %(debug purposes)
    % Initial data:
    h = a*exp(-(Zeta-4).^2/ (2* kappa^4 * sigma^2));
    u = -h; % necessary velocity for unidirectional solution (down-going mode only)
    v = zeros(size(h));

    %J = ones(size(J));
end

h = neumann_correction(h);
u = impermiability_u(u);
v = impermiability_v(v);

[h, u, v] = enforce_barrier(h, u, v, th_zeta, th_xi);  % Enforce slit barrier

%% Setting wave data

H = reshape(h, numel(data.J),1);
U = reshape(u, numel(data.J),1);
V = reshape(v, numel(data.J),1);

%% Precompute the forward mapping
%z = eval(f, w);

%{
%z1 = eval(f, w1);
%z2 = eval(f, w2);
%z3 = eval(f, w3);
%}

z1=z(:,1:th_xi); %branch 1 
z2=z(1:th_zeta,th_xi:end); %branch 2
z3=z(th_zeta+1:end,th_xi:end); %branch 3

% Extract the real and imaginary parts of the forward-mapped grid
X = real(z);
Y = imag(z);

X1=real(z1); Y1=imag(z1);
X2=real(z2); Y2=imag(z2);
X3=real(z3); Y3=imag(z3);

%% Plot initial data

surf(Xi, Zeta, h), hold on
plot3([data.xi(th_xi) data.xi(end)], [data.zeta(th_zeta) data.zeta(th_zeta)], [0 0], 'r-', 'LineWidth', 2);


%sz = size(J);
%J = ones(sz);

%% Main loop % RK4 time-stepping
%dist = 0;
t = 0; iter=0;
t_array = 0:dt:T;
tol = 0.00001*a;
%while t < T
%while dist < options.travel_distance*comp_efetivo_can
while max(h(:, end)) < tol
    %h_pre = h;
    %u_pre = u;
    %v_pre = v;
    
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
    
    v = impermiability_v(v);
    h = neumann_correction(h);
    u = impermiability_u(u);
    
    [h, u, v] = enforce_barrier(h, u, v, th_zeta, th_xi);  % Enforce slit barrier
    
    surf(Xi, Zeta, h)
    
    %% Saves selected number of frames
    if mod(iter,floor(numel(t_array)/options.frames))==0 || t == T-1-dt
        
        H = [H, reshape(h,numel(J),1)];
        U = [U, reshape(u,numel(J),1)];
        V = [V, reshape(v,numel(J),1)];
        
        surf(Xi, Zeta, h), hold on
        plot3([data.xi(th_xi) data.xi(end)], [data.zeta(th_zeta) data.zeta(th_zeta)], [0 0], 'r-', 'LineWidth', 2);
        
        
        %plot3(data.xi(th_xi+10),data.zeta(th_zeta), [0 0], 'x')
        %plot3(data.xi(th_xi+10),data.zeta(th_zeta+1), [0 0], 'x')
        %plot3(data.xi(th_xi+10),data.zeta(th_zeta-1), [0 0], 'x')
        %plot3(data.xi(th_xi+10),data.zeta(th_zeta-2), [0 0], 'x')
        
        hold off
        
    end
    
    %{
    if mod(iter,200)==0 && strcmp(visual,'physical')
    
        %saveframe()
         
    hh=h;
    
    jmp = 1;
    
    h1=hh(1:end,1:jmp:th_xi+1);
    h2=hh(1:th_zeta,th_xi-1:jmp:end);
    h3=hh(th_zeta+1:end,th_xi-1:jmp:end);
    
    XX1 = X1(1:end,1:jmp:end);
    YY1 = Y1(1:end,1:jmp:end);
    
    XX2 = X2(1:end,1:jmp:end);
    YY2 = Y2(1:end,1:jmp:end);
    
    XX3 = X3(1:jmp:end,1:jmp:end);
    YY3 = Y3(1:jmp:end,1:jmp:end);
    
    %subplot(1,2,1)
    mesh(XX1, YY1, h1, 'edgecolor', 'k'); hold on,
    mesh(XX2, YY2, h2, 'edgecolor', 'k');
    mesh(XX3, YY3, h3, 'edgecolor', 'k');
    hold off,
    view(az, el);
    zlim([-0.15,.2])
    %caxis([min(h(:)), max(h(:))]);  % Set the color axis limits based on the data range
    xlabel('X'); ylabel('Y'); zlabel('h');
    title(['Time evolution of wave profile = ',num2str(t)]);
    
    %subplot(1,2,2)
    %surf(Xi,Zeta,h)
    
    drawnow;
    
    elseif mod(iter,50)==0 && strcmp(visual,'canonical')
        %
        hh=h;
        
        h1=hh(:,1:th_xi+1);
        h2=hh(1:th_zeta,th_xi-1:end);
        h3=hh(th_zeta+1:end,th_xi-1:end);
        
        %subplot(1, 2, 1);
        %mesh(X1, Y1, h1, 'edgecolor', 'k'); hold on,
        %mesh(X2, Y2, h2, 'edgecolor', 'k');
        %mesh(X3, Y3, h3, 'edgecolor', 'k');
        
        %hold off,
        %
        %subplot(1, 2, 2);
        mesh(xi,zeta,h)
        zlim([-0.05,a])
        
        %set(gcf, 'Renderer', 'opengl');  % Better rendering quality
        %set(gca, 'FontSize', 14);        % Make axes labels crisper
        %shading interp                   % Smooth surface shading
        %lighting gouraud                 % Optional: if you use lighting
        
        drawnow;
        frame = getframe(gcf); % Capture current figure
        writeVideo(vwriter, frame); % Write frame to video
        
    end
  
    pause(0.001)
    %}
    % Update time
    t = t + dt;
    %hh = h(floor(end/2),:);   
    %[~,loc] = findpeaks(hh,'MinPeakHeight', 0.7*a); % finds index for which final wave prof peaks (center of final gaussian)
    %x0f = data.xi(loc);
    %dist = abs(xi0 - x0f); %gets distance between current and initial peak positions
    
    iter=iter+1;
    max(h(:, end))
    
    if t > T
        message = ['Exceeded final time ceiling. T = ', num2str(T)]
        break
    end
end

%% Save data if requested
    if options.want_save
        ang_display = round(angles, 3);
        save(['WaveData/kappa', num2str(kappa),'widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat'])
    end
  
% Aux functions    
   
%% Slit boundary condition function
function [h,u, v] = enforce_barrier(h, u, v, b1, b2)
    % Strict barrier enforcement
    h(b1+1,b2:end) = h(b1+2,b2:end);   % h continues smoothly
    u(b1+1,b2:end) = 0;
    %v(b1+1,b2:end) = 0;
    
    h(b1-1, b2:end) = h(b1, b2:end); % h continues smoothly
    u(b1,b2:end) = 0;
    %v(b1,b2:end) = 0;
end


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

