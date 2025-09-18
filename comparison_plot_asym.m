% Create figure with different wave data, symmetric case
clear all, close all, clc

addpath('External/')

%Setting parameters
widths = [5, 5, 5]; 

parameter_station % Go through preferred secondary arguments

kappas = [0.2,0.3];
thetas3 = 5*pi/12;

%% Parameter sweep color grid view
anglesStart = [0, pi - pi/30, pi + 0];
angles = anglesStart;


D = zeros(length(kappas),length(thetas3));

for ii = 1:length(kappas)
    kappa = kappas(ii);
    lambda_f = widths(1)/kappa;
    Lx = lambda_f * (travel_distance + 1) / 2;
    
    for jj = 1:length(thetas3)
        theta = thetas3(jj);
        angles(3) = anglesStart(3) + theta;
    
        [h1,h2,h3,th_xi,xi] = loading(kappa, angles, widths,Lx);

        D(ii,jj) = (max(h3)-max(h2))/abs(max(h3));

        if deltaHeightPlot
            Dindex = (ii-1)*length(thetas3) + jj;

            xi3 = xi(th_xi:end);

            tmp = xi3(1);
            %dilation = (s/l);
            dilation = 1;

            xi3 = xi3 - tmp;
            xi3 = dilation*xi3;
            xi3 = xi3 + tmp;

            % call helper for each branch
            [loc2, pk2] = find_main_peak(h2);
            [loc3, pk3] = find_main_peak(h3);
            
            xpeak2 = xi3(loc2);
            xpeak3 = xi3(loc3);
            shift = xpeak2 - xpeak3;   % amount to move branch k so peaks coincide in x

            xi3_aligned = xi3 + shift;
            
            subplot(length(kappas),length(thetas3),Dindex) 
            plot(xi3_aligned, h3, '-','LineWidth',2.2,...
                'Color', 'b', 'DisplayName','Branch 2'); hold on
            plot(xi3, h2, '-.','LineWidth',2.2,...
                'Color', 'r', 'DisplayName','Branch 3'); hold off
            
            title(['kappa= ', num2str(kappa),' theta= ', num2str(theta)])
            legend('show')

            xlabel('\xi')

            set(gca, 'FontSize',18, ...      % tick labels larger
                 'LineWidth',1.5, ...    % axis lines thicker
                 'TickDir','out', ...    % ticks outward
                 'Box','off')            % remove top/right frame
            

        end
    
    end
    angles = anglesStart;

end


if parameterSweep
figure(2)
%imagesc(D)
imagesc(thetas3,kappas,D)
ylabel('\kappa', 'Rotation',0)
xlabel('\theta_{asym}')
colorbar

set(gca, 'FontSize',18, ...      % tick labels larger
     'LineWidth',1.5, ...    % axis lines thicker
     'TickDir','out', ...    % ticks outward
     'Box','off')            % remove top/right frame
end


% Aux functions
%{

%% Transmitted waves comparison
figure(1)
for ii = 1:length(angles)
    for jj = 1:length(kappas)
        [~,h2,h3, th_xi,xi] = loading(kappas{jj},angles{ii}, widths);
        
    end
end

% low kappa, low angle
[~,h2,h3, th_xi,xi] = loading(kappa_low,angles_low, widths);
subplot(3,2,1) 
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'Color', 'b', 'DisplayName', ['Angle = ',num2str(angles_low(2))]), hold on 
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'Color', 'r','DisplayName', ['Angle = ',num2str(angles_low(3))]), hold off
title('Kappa low, angle low')

% low kappa, high angle
[~,h2,h3, th_xi,xi] = loading(angles_low, kappa_high, widths);
subplot(3,2,5) 
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'Color', 'b', 'DisplayName', ['Angle = ',num2str(angles_low(2))]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'Color', 'r','DisplayName', ['Angle = ',num2str(angles_low(3))]), hold off
title('Kappa high, angle low')

[~,h2,h3, th_xi,xi] = loading(angles_low, kappa_mid, widths);
subplot(3,2,3) % mid kappa, low angle
plot(xi(th_xi:end), h2, '-.','LineWidth',2,'Color', 'b', 'DisplayName', ['Angle = ',num2str(angles_high(2))]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2,'Color', 'r', 'DisplayName', ['Angle = ',num2str(angles_high(3))]), hold off
title('Kappa low, angle high')

[~,h2,h3, th_xi,xi] = loading(angles_low, kappa_high, widths);
subplot(3,2,4) % high kappa, low angle
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'Color', 'b', 'DisplayName', ['Angle = ',num2str(angles_low(2))]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'Color', 'r','DisplayName', ['Angle = ',num2str(angles_low(3))]), hold off
title('Kappa high, angle low')

[~,h2,h3, th_xi,xi] = loading(angles_high, kappa_low, widths);
subplot(3,2,2) % low kappa, high angle
plot(xi(th_xi:end), h2, '-.','LineWidth',2,'Color', 'b', 'DisplayName', ['Angle = ',num2str(angles_high(2))]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2,'Color', 'r', 'DisplayName', ['Angle = ',num2str(angles_high(3))]), hold off
title('Kappa low, angle high')

[~,h2,h3, th_xi,xi] = loading(angles_high, kappa_high, widths);
subplot(3,2,6) % high kappa, high angle
plot(xi(th_xi:end), h2, '-.','LineWidth',2,'Color', 'b', 'DisplayName', ['Angle = ',num2str(angles_high(2))]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'Color', 'r','DisplayName', ['Angle = ',num2str(angles_high(3))]), hold off
title('Kappa high, angle high')
%}
%---------------------
function [h1,h2,h3,th_xi,xi] = loading(kappa, angles, widths, Lx)
ang_display = round(angles .* 1000) ./ 1000;
load(['GraphData/widths=', mat2str(widths), '_angles=', mat2str(ang_display),'_length=',...
            mat2str(Lx), '.mat'],...
        'w','th_xi','th_zeta');

load(['WaveData/kappa=', num2str(kappa),'_widths=', mat2str(widths), '_angles=', mat2str(ang_display), '_length=',...
            mat2str(Lx),'.mat'],'h')

%tmp = size(H);

%h = reshape(H(:,tmp(2)),size(z));
h1 = h(floor(end/2),1:th_xi);
h2 = h(floor((end+th_zeta)/2),th_xi:end);
h3 = h(floor((1+th_zeta)/2),th_xi:end);

xi = real(w);
xi = xi(1,:);

end

function [loc,pkval] = find_main_peak(h) %GPT suggested function (overkill, could be simpler)
    % try plain findpeaks (largest)
    [pks,locs] = findpeaks(h, 'SortStr','descend');
    if ~isempty(locs)
        loc = locs(1);
        pkval = pks(1);
        return
    end
    % try with a modest prominence threshold
    prom = 0.1*(max(h)-min(h));
    [pks,locs] = findpeaks(h, 'MinPeakProminence', prom, 'SortStr','descend');
    if ~isempty(locs)
        loc = locs(1);
        pkval = pks(1);
        return
    end
    % fallback: global maximum
    [pkval, loc] = max(h);
end


