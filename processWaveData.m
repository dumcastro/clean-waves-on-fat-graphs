% Copyright (C) 2025 Eduardo Castro
% Author: Eduardo Castro <eduardocastro@MacBook-Air-de-Eduardo.local>
% Created: 2025-06-03

function [] = processWaveData(kappa, widths, angles, options)


% Set default parameter values
    if nargin < 4
        options = struct();
    end

    % Default graph parameters
    if nargin < 1 || isempty(kappa), kappa = 0.25; end
    if nargin < 2 || isempty(widths), widths = [1, 0.5, 0.5]; end
    if nargin < 3 || isempty(angles), angles = [0, 2*pi/3, 4*pi/3]; end

    % Default numerical parameters
    default_options = struct(...
        'plot_physical', false,...
        'plot_canonical', false,...
        'play_movie_phys', false,...
        'play_movie_canonical', true,...
        'twoD_plot', false,...
        'twoD_animation', true,...
        'save_video', false,...
        'jmp_xi', 1,...
        'jmp_zeta', 1,...
        'az',-20,...
        'el',40);
    

    % Merge user options with defaults
    option_names = fieldnames(default_options);
    for k = 1:length(option_names)
        if ~isfield(options, option_names{k})
            options.(option_names{k}) = default_options.(option_names{k});
        end
    end

    %% Load wave Graph
    ang_display = round(angles .* 1000) ./ 1000;
    data = load(['WaveData/kappa', num2str(kappa),'widths= ',...
        mat2str(widths), 'angles= ', mat2str(ang_display), '.mat']);

    %%
    h = data.h;
    %t = data.t;
    z = data.z;
    %xi_lims = data.xi_lims;
    %Xi = data.Xi;
    %J = data.J;
    th_zeta = data.th_zeta;
    th_xi = data.th_xi;

    
    %% Dividing the domain into sectors from plotting purposes
    
    jmp = options.jmp_xi;
    jmpz = options.jmp_zeta;

    ngrid = size(h);

    [Nzeta, ~] = deal(ngrid(1),ngrid(2));
    
    if numel(widths) == 3
    zindexes1 = [2:jmpz:th_zeta-1, th_zeta, th_zeta+1, th_zeta+2:jmpz:Nzeta-1];
    zindexes2 = [2:jmpz:th_zeta-1, th_zeta];
    zindexes3 = [th_zeta+1, th_zeta+2:jmpz:Nzeta-1];
    
    z1=z(zindexes1,[1:jmp:th_xi,th_xi]); %branch 1 
    z2=z(zindexes2,th_xi:jmp:end); %branch 2
    z3=z(zindexes3,th_xi:jmp:end); %branch 3

    % Extract the real and imaginary parts of the forward-mapped grid
    %X = real(z);
    %Y = imag(z);

    X1=real(z1); Y1=imag(z1);
    X2=real(z2); Y2=imag(z2);
    X3=real(z3); Y3=imag(z3);
    end
    
    %%
    %{
    if options.plot_physical %(This is outdated)
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
    view(options.az, options.el);
    zlim([-0.15,.2])
    %caxis([min(h(:)), max(h(:))]);  % Set the color axis limits based on the data range
    xlabel('X'); ylabel('Y'); zlabel('h');
    title(['Time evolution of wave profile = ',num2str(t)]);

    %subplot(1,2,2)
    %surf(Xi,Zeta,h)

    drawnow;
    end
    %}


    if options.plot_canonical
        %

        hh=h;

        h1=hh(:,1:th_xi);
        h2=hh(1:th_zeta,th_xi:end);
        h3=hh(th_zeta+1:end,th_xi:end);

        %subplot(1, 2, 1);
        %mesh(X1, Y1, h1, 'edgecolor', 'k'); hold on,
        %mesh(X2, Y2, h2, 'edgecolor', 'k');
        %mesh(X3, Y3, h3, 'edgecolor', 'k');

        %hold off,
        %
        %subplot(1, 2, 2);
        
        
        
        mesh(real(data.data.w),imag(data.data.w),h)
        zlim([-0.05,a])

        %set(gcf, 'Renderer', 'opengl');  % Better rendering quality
        %set(gca, 'FontSize', 14);        % Make axes labels crisper
        %shading interp                   % Smooth surface shading
        %lighting gouraud                 % Optional: if you use lighting

        drawnow;

    end
    
    tmp = size(data.H);
    
    %% Play movie and save video
    if options.play_movie_phys
        if options.save_video
            vwriter = VideoWriter(['Export/kappa',...
                num2str(kappa),'widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mp4'], 'MPEG-4');
            vwriter.FrameRate = 10;      % Adjust for desired smoothness
            vwriter.Quality = 100;       % Max quality (optional)
            open(vwriter);
        end
        

        %mytitle = ['Angle = ', num2str(rad2deg(angles(3)-angles(2))), ' degrees'];
        
        
        figure
        %for i = 1:data.options.frames-2
        for i = 1:tmp(2)

            h = reshape(data.H(:,i),size(z));

            if numel(widths) == 3

            hh=h;

            %{
            jmp = options.jmp_xi;
            jmpz = options.jmp_zeta;
            
            ngrid = size(h);
            
            [Nzeta, ~] = deal(ngrid(1),ngrid(2));
            
            zindexes1 = [2:jmpz:th_zeta-1, th_zeta, th_zeta+1, th_zeta+2:jmpz:Nzeta-1];
            zindexes2 = [2:jmpz:th_zeta-1, th_zeta];
            zindexes3 = [th_zeta+1, th_zeta+2:jmpz:Nzeta-1];
            %}
            h1=hh(zindexes1,[1:jmp:th_xi,th_xi]);
            h2=hh(zindexes2,th_xi:jmp:end);
            h3=hh(zindexes3,th_xi:jmp:end);
            
            %{
            XX1 = X1(zindexes1,[1:jmp:end, end]);
            YY1 = Y1(zindexes1,[1:jmp:end, end]);
            
            XX2 = X2(zindexes2,1:jmp:end);
            YY2 = Y2(zindexes2,1:jmp:end);

            XX3 = X3(zindexes3,1:jmp:end);
            YY3 = Y3(zindexes3,1:jmp:end);
            %}
            %% Debug
            %---------
            %h1 = zeros(size(h1));
            %h2 = zeros(size(h2));
            %---------
            mesh(X1, Y1, h1, 'edgecolor', 'k'); hold on,
            mesh(X2, Y2, h2, 'edgecolor', 'k');
            mesh(X3, Y3, h3, 'edgecolor', 'k');   
            
            %mesh(data.X(2:end-1,:),data.Y(2:end-1,:),h(2:end-1,:))
            
            hold off,
            view(options.az, options.el);
            zlim([-0.02,.12])
            %caxis([min(h(:)), max(h(:))]);  % Set the color axis limits based on the data range
            xlabel('X'); ylabel('Y'); zlabel('h','Rotation', 0);
            %title(['Time evolution of wave profile = ',num2str(t)]);
            %title(mytitle)

            set(gca, 'FontSize', 16)   % makes axis numbers larger
      
            pause(0.3)

            drawnow;
            else
                mesh(data.X,data.Y,h)
                view(options.az, options.el);
                %zlim([-0.02,.12])
                xlabel('X'); ylabel('Y'); zlabel('h','Rotation', 0);
                %title(['Time evolution of wave profile = ',num2str(t)]);
                %title(mytitle)
    
                set(gca, 'FontSize', 16)   % makes axis numbers larger
          
                pause(0.2)
    
                drawnow;
            end

            
            if options.save_video
                %while counter/i < 10
                frame = getframe(gcf); % Capture current figure
                writeVideo(vwriter, frame); % Write frame to video
                %counter = counter + 1;
                %end
            end

        end

    end
    
    
    %% Animation in canonical coordinates
    if options.play_movie_canonical
        
        %mytitle = ['Angle = ', num2str(rad2deg(angles(3)-angles(2))), ' degrees: Canonical domain'];
   
        figure

        %for i = 1:data.options.frames-2
        for i = 1:tmp(2)
        h = reshape(data.H(:,i),size(z));
            
        %hh=h;

        %h1=hh(:,1:th_xi+1);
        %h2=hh(1:th_zeta,th_xi-1:end);
        %h3=hh(th_zeta+1:end,th_xi-1:end);

        %subplot(1, 2, 1);
        %mesh(X1, Y1, h1, 'edgecolor', 'k'); hold on,
        %mesh(X2, Y2, h2, 'edgecolor', 'k');
        %mesh(X3, Y3, h3, 'edgecolor', 'k');

        %hold off,
        %
        %subplot(1, 2, 2);
        mesh(real(data.data.w),imag(data.data.w),h)
        %zlim([-0.05,0.1])

        
        %title(mytitle)
        drawnow;
        

        pause(0.2)
        end
        
    end
    

    %% 2D plot
    if options.twoD_plot
        if numel(widths) == 2
        
            h = reshape(data.H(:,tmp(2)),size(z));
            h = h(floor(end/2),:);   % final snapshot
            
            hh = reshape(data.H(:,1),size(z));
            hh = hh(floor(end/2),:); % initial snapshot
            
            xi = real(data.data.w);   % x-axis
            xi = xi(1,:);
            
            % Find main peaks
            [~,loc_h]  = findpeaks(h, xi, 'SortStr','descend', 'NPeaks',1);
            [~,loc_hh] = findpeaks(hh, xi, 'SortStr','descend', 'NPeaks',1);

            % Compute shift needed
            shift = loc_h - loc_hh;
            
            % Shift the initial profile
            hh_shifted = interp1(xi, hh, xi - shift, 'linear', 0);

            figure;

            % Use strong, distinguishable colors
            plot(xi, h, 'LineWidth', 2.5, 'Color',[0 0.45 0.74])        % MATLAB blue
            hold on
            plot(xi, hh_shifted, '--', 'LineWidth', 2.5, 'Color',[0.85 0.33 0.1]) % MATLAB orange
            
            legend({'Final snapshot','Initial snapshot shifted'}, ...
                   'Location','best', 'FontSize',18, 'Box','off')
            
            xlabel('\xi','FontSize',20)
            %ylabel('Amplitude','FontSize',20)
            
            set(gca, 'FontSize',18, ...      % tick labels larger
                     'LineWidth',1.5, ...    % axis lines thicker
                     'TickDir','out', ...    % ticks outward
                     'Box','off')            % remove top/right frame



        elseif numel(widths) == 3
            
            h = reshape(data.H(:,tmp(2)),size(z));
            h1 = h(floor(end/2),1:th_xi);
            h2 = h(floor((1+th_zeta)/2),th_xi:end);
            h3 = h(floor((end+th_zeta)/2),th_xi:end);
            
            
            xi = real(data.data.w);
            xi1 = xi(floor(end/2),1:th_xi);
            xi2 = xi(floor((1+th_zeta)/2),th_xi:end);
            xi3 = xi(floor((end+th_zeta)/2),th_xi:end);
            
            
            figure;
            subplot(3, 1, 1);
            plot(xi1,h1)
            %title('Physical Region');

            subplot(3, 1, 2);
            plot(xi2,h2)
            %title('Numerical canonical domain');

            subplot(3, 1, 3);
            plot(xi3,h3)
            %title('Canonical domain width = 1');
   
        end   
    end
    
    %% 2D animation on each branch
    if options.twoD_animation
       
        %mytitle = ['Angle = ', num2str(rad2deg(angles(3)-angles(2))), ' degrees'];
        
        
        figure
        %for i = 1:data.options.frames-2
        for i = 1:tmp(2)

            h = reshape(data.H(:,i),size(z));

            %{
            jmp = options.jmp_xi;
            jmpz = options.jmp_zeta;
            
            ngrid = size(h);
            
            [Nzeta, ~] = deal(ngrid(1),ngrid(2));
            
            zindexes1 = [2:jmpz:th_zeta-1, th_zeta, th_zeta+1, th_zeta+2:jmpz:Nzeta-1];
            zindexes2 = [2:jmpz:th_zeta-1, th_zeta];
            zindexes3 = [th_zeta+1, th_zeta+2:jmpz:Nzeta-1];
            %}
            
            h1 = h(floor(end/2),1:data.th_xi);
            h2 = h(floor((1+th_zeta)/2),th_xi:end);
            h3 = h(floor((end+th_zeta)/2),th_xi:end);
            
            
            xi = real(data.data.w);
            xi1 = xi(floor(end/2),1:th_xi);
            xi2 = xi(floor((1+th_zeta)/2),th_xi:end);
            xi3 = xi(floor((end+th_zeta)/2),th_xi:end);
            

            subplot(2,1,1) % branch i
            plot(xi1, h1,'LineWidth',1,'DisplayName','Branch i')
            title('Reflected wave') 
            legend('show')
            
            subplot(2,1,2) % branch j
            plot(xi2, h2, '-.','LineWidth',1,'DisplayName','Branch j'), hold on
            plot(xi3, h3, '--','LineWidth',1,'DisplayName','Branch k'), hold off
            title('Transmitted waves') 
            legend('show')
            
            pause(0.3)

            drawnow;
            
        end

    end
    
    if options.save_video
       close(vwriter)
    end


