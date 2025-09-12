% Create figure with different wave data, symmetric case
clear all, close all, clc

addpath('External/')

%Setting parameters
widths = [5, 5, 5]; 

angles_high = [0, pi - pi/24, pi + pi/2 - pi/12];
angles_low = [0, pi - pi/24, pi + pi/12];

theta_high = num2str(rad2deg(angles_high(3)-angles_high(2)));
theta_low = num2str(rad2deg(angles_low(3)-angles_low(2)));

kappa_high = 0.3;
kappa_low = 0.1;

parameter_station % Go through preferred secondary arguments

%angles = angles_low;
%kappa = kappa_low;

figure(1)

%% Transmitted waves comparison
%subplot(2,2,1) % low kappa, low angle
[~,h2,h3, th_xi,xi] = loading(angles_low, kappa_low, widths);
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'DisplayName', ['Angle = ',angles_low(2)]), hold on 
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'DisplayName', ['Angle = ',angles_low(3)]), hold off

%processWaveData(kappa_low, widths, angles_low,wave_vis_options)

%{
subplot(2,2,2) % high kappa, low angle
[~,h2,h3, th_xi,xi] = loading(angles_low, kappa_high, widths);
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'DisplayName', ['Angle = ',angles_low(2)]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'DisplayName', ['Angle = ',angles_low(3)]), hold off

subplot(2,2,3) % low kappa, high angle
[~,h2,h3, th_xi,xi] = loading(angles_high, kappa_low, widths);
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'DisplayName', ['Angle = ',angles_high(2)]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'DisplayName', ['Angle = ',angles_high(3)]), hold off

subplot(2,2,4) % high kappa, high angle
[~,h2,h3, th_xi,xi] = loading(angles_high, kappa_high, widths);
plot(xi(th_xi:end), h2, '-.','LineWidth',2, 'DisplayName', ['Angle = ',angles_high(2)]), hold on
plot(xi(th_xi:end), h3, '--','LineWidth',2, 'DisplayName', ['Angle = ',angles_high(3)]), hold off
%}

%---------------------
function [h1,h2,h3,th_xi,xi] = loading(angles, kappa, widths)
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


