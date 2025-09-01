% Create figure with different wave data, symmetric case
clear all, close all, clc

addpath('External/')

%Setting parameters
widths = [2.5, 2.5]; %fix

angles_high = [0, pi + pi/2 - pi/24];
angles_low = [0, pi + pi/24];

theta_high = '90 degrees';
theta_low = '7.5 degrees';

kappa_high = 0.4;
kappa_mid = 0.1;
kappa_low = 0.04;

parameter_station % Go through preferred secondary arguments

%angles = angles_low;
%kappa = kappa_low;


%Go for it...!

figure(1)

% Low kappa
subplot(2,1,1) % low kappa
[h1,~,xi] = loading(angles_low, kappa_low, widths);
plot(xi, h1,'LineWidth',2.5, 'DisplayName', ['Angle = ',theta_low]), hold on
%xlim([0,180])

[h1,~,xi] = loading(angles_high, kappa_low, widths);
plot(xi, h1, '-.','LineWidth',2.5, 'DisplayName', ['Angle = ',theta_high]), hold off
legend('show') % Show legend
title(['\kappa = ', num2str(kappa_low)]) % Title for Low Kappa subplot

set(gca, 'FontSize',18, ...      % tick labels larger
                     'LineWidth',2.5, ...    % axis lines thicker
                     'TickDir','out', ...    % ticks outward
                     'Box','off')            % remove top/right frame

% Mid kappa
subplot(2,1,2) % mid kappa
[h1,~,xi] = loading(angles_low, kappa_mid, widths);
plot(xi, h1,'LineWidth',2.5, 'DisplayName', ['Angle = ',theta_low]), hold on

[h1,~,xi] = loading(angles_high, kappa_mid, widths);
plot(xi, h1, '-.','LineWidth',2.5,'DisplayName', ['Angle = ',theta_high]), hold off
%legend('show') % Show legend
title(['\kappa = ', num2str(kappa_mid)]) % Title for Mid Kappa subplot

set(gca, 'FontSize',18, ...      % tick labels larger
                     'LineWidth',1.5, ...    % axis lines thicker
                     'TickDir','out', ...    % ticks outward
                     'Box','off')            % remove top/right frame



%{
% High kappa
subplot(3,1,3) % high kappa
[h1,~,xi] = loading(angles_low, kappa_high, widths);
plot(xi, h1, 'LineWidth',1,'DisplayName', ['Angle = ',theta_low]), hold on
%xlim([0,30])

[h1,~,xi] = loading(angles_high, kappa_high, widths);
plot(xi, h1, '-.','LineWidth',1,'DisplayName', ['Angle = ',theta_high]), hold off
%legend('show') % Show legend
title(['\kappa = ', num2str(kappa_high)]) % Title for High Kappa subplot
%}

%---------------------
function [h1,th_xi,xi] = loading(angles, kappa, widths)

%lambda_f = widths(1)/kappa; %fix
%travel_distance = 2; %fix
%Lx = lambda_f * (travel_distance + 1) / 2; %fix

ang_display = round(angles .* 1000) ./ 1000;
data = load(['WaveData/kappa', num2str(kappa),'widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat']);

th_xi = data.th_xi;
tmp = size(data.H);

h = reshape(data.H(:,tmp(2)),size(data.z));
h1 = h(floor(end/3),:);
%h1 = h1(1:data.th_xi);

xi = real(data.data.w);
xi = xi(1,:);


end


