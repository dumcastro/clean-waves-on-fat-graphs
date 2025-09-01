clear all, clc, close all

widths = [5, 5];
widths2 = 2*[5, 5];
angles = [0, pi+pi/6];
kappa = 0.1;
parameter_station

ang_display = round(angles .* 1000) ./ 1000;

data = load(['WaveData/kappa', num2str(kappa),...
    'widths= ', mat2str(widths), 'angles= ', mat2str(ang_display), '.mat']);

data2 = load(['WaveData/kappa', num2str(kappa),...
    'widths= ', mat2str(widths2), 'angles= ', mat2str(ang_display), '.mat']);

%h = data.h;
z = data.z;
z2 = data2.z;



tmp = size(data.H);
tmp2 = size(data2.H);


h = reshape(data.H(:,tmp(2)),size(z));
h = h(floor(end/2),:);   % final snapshot

h2 = reshape(data2.H(:,tmp2(2)),size(z2));
h2 = h2(floor(end/2),:);   % final snapshot

%h2 = h2(1:2:end-30); %shaving ?


xi = real(data.data.w);   % x-axis
xi = xi(1,:);

xi2 = real(data2.data.w);   % x-axis
xi2 = xi2(1,:);

plot(2*xi, h, 'LineWidth', 2.5, 'Color',[0 0.45 0.74]), hold on
plot(xi2, h2, '--', 'LineWidth', 2.5, 'Color',[0.85 0.33 0.1])

legend({'\lambda = 50, w = 5','\lambda=100, w = 10'}, ...
                   'Location','best', 'FontSize',18, 'Box','off')

set(gca, 'FontSize',18, ...      % tick labels larger
                     'LineWidth',1.5, ...    % axis lines thicker
                     'TickDir','out', ...    % ticks outward
                     'Box','off')            % remove top/right frame



