% Create figure with different wave data, symmetric case
clear all, close all, clc

addpath('External/')

%Setting parameters
widths = [5, 5, 5]; 

parameter_station % Go through preferred secondary arguments

anglesStart = [0, pi - pi/30, pi + 0];
angles = anglesStart;

prof_plot = false;

%{
angles_high = [0, pi - pi/24, pi + pi/2 - pi/12];
angles_low = [0, pi - pi/24, pi + pi/12];

theta_high = num2str(rad2deg(angles_high(3)-angles_high(2)));
theta_low = num2str(rad2deg(angles_low(3)-angles_low(2)));

kappa_high = 0.4;
kappa_mid = 0.25;
kappa_low = 0.1;

angles = {angles_low, angles_high};
kappas = {kappa_low, kappa_mid, kappa_high};

%angles = angles_low;
%kappa = kappa_low;
%}

kappas = 0.1:0.1:0.4;
thetas = pi/12: pi/6 :pi/2-pi/12;

D = zeros(length(kappas),length(thetas));

for ii = 1:length(kappas)
    kappa = kappas(ii);
    
    for jj = 1:length(thetas)
        theta = thetas(jj);
        angles(3) = anglesStart(3) + theta;
    
        [h1,h2,h3,th_xi,xi] = loading(kappa, angles, widths);

        D(ii,jj) = (max(h3)-max(h2))/abs(max(h3));


        if prof_plot
            Dindex = (ii-1)*length(thetas) + jj;

            if Dindex == 6
                figure

                xi3 = xi(th_xi:end);

                tmp = xi3(1);
                %dilation = (s/l);
                dilation = 0.645;
    
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

                plot(xi3_aligned, h3, '-','LineWidth',2.2,...
                    'Color', 'b', 'DisplayName','Branch k (shifted)'); hold on
                plot(xi3, h2, '-.','LineWidth',2.2,...
                    'Color', 'r', 'DisplayName','Branch j'); hold off
                
                title('Transmitted waves (aligned peaks)')
                legend('show')

                xlabel('\xi')

                set(gca, 'FontSize',18, ...      % tick labels larger
                     'LineWidth',1.5, ...    % axis lines thicker
                     'TickDir','out', ...    % ticks outward
                     'Box','off')            % remove top/right frame

                p = 32;
            end

            subplot(length(kappas),length(thetas),Dindex) 
            plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'Color', 'b'), hold on 
            plot(xi(th_xi:end), h3, '--','LineWidth',2, 'Color', 'r'), hold off
            title(['kappa= ', num2str(kappa),' theta= ', num2str(theta)])

        end
    
    end
    angles = anglesStart;

end

figure(2)
%imagesc(D)
imagesc(thetas,kappas,D)
ylabel('\kappa', 'Rotation',0)
xlabel('\theta_{asym}')
colorbar

set(gca, 'FontSize',18, ...      % tick labels larger
     'LineWidth',1.5, ...    % axis lines thicker
     'TickDir','out', ...    % ticks outward
     'Box','off')            % remove top/right frame


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
function [h1,h2,h3,th_xi,xi] = loading(kappa, angles, widths)
ang_display = round(angles .* 1000) ./ 1000;
data = load(['WaveData/kappa', num2str(kappa),'widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat']);

th_xi = data.th_xi;
th_zeta = data.th_zeta;
tmp = size(data.H);

h = reshape(data.H(:,tmp(2)),size(data.z));
h1 = h(floor(end/2),1:data.th_xi);
h2 = h(floor((end+th_zeta)/2),th_xi:end);
h3 = h(floor((1+th_zeta)/2),th_xi:end);

xi = real(data.data.w);
xi = xi(1,:);


end

function [loc,pkval] = find_main_peak(h)
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


