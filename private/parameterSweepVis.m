function [] = parameterSweepVis(widthTriplets, kappas, theta2, thetas3,travelDistance,options)

anglesStart = [0, pi - theta2, pi + 0];
angles = anglesStart;

DFixedWidths = zeros(length(kappas),length(thetas3));
DFixedTheta = zeros(length(kappas),length(widthTriplets));
DFixedKappa = zeros(length(widthTriplets),length(thetas3));

widthsAsymRatio = [];

for kk = 1:length(widthTriplets)
    widths = widthTriplets{kk};
    widthsAsymRatio = [widthsAsymRatio, widths(3)/widths(2)];

for ii = 1:length(kappas)
    kappa = kappas(ii);
    lambda_f = widths(1)/kappa;
    Lx = lambda_f * (travelDistance + 1) / 2;
    
    for jj = 1:length(thetas3)
        theta = thetas3(jj);
        angles(3) = anglesStart(3) + theta;
    
        [~,h2,h3,th_xi,xi] = loading(kappa, angles, widths,Lx);

        DFixedWidths(ii,jj) = (max(h3)-max(h2))/abs(max(h3));
        DFixedTheta(ii,kk) = (max(h3)-max(h2))/abs(max(h3));
        DFixedKappa(kk,jj) = (max(h3)-max(h2))/abs(max(h3));


        if options.deltaHeightPlot % (for fixed widths only)
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
end

%% theta and kappa influence on deltaHeight 
if options.colorGridWidthsFixed

    if ~(all(widthTriplets{1} == [5,5,5]) && length(widthTriplets)==1)
        error('Adjust width to 5 5 5')
    end

    figure(2)
    %imagesc(D)
    imagesc(thetas3,kappas,DFixedWidths)
    ylabel('\kappa', 'Rotation',0)
    
    yticks(kappas)
    
    xlabel('\theta_{asym}')
    
    xticks(thetas3)
    
    xticklabels({'\pi/12','3\pi/12','5\pi/12'})
    
    colorbar
    
    set(gca, 'FontSize',18, ...      % tick labels larger
         'LineWidth',1.5, ...    % axis lines thicker
         'TickDir','out', ...    % ticks outward
         'Box','off')            % remove top/right frame
end

%% theta and width influence...

%
if options.colorGridKappaFixed

    if ~ isscalar(kappas)
        error('Choose a fixed kappa')
    else
        disp(['kappa is ', num2str(kappa)])
    end

    figure(3)
    %imagesc(D)
    imagesc(thetas3,widthsAsymRatio,DFixedKappa)
    ylabel('widths', 'Rotation',0)
    
    yticks(widthsAsymRatio)
    
    xlabel('\theta_{asym}')
    
    xticks(thetas3)
    
    %xticklabels({'\pi/12','3\pi/12','5\pi/12'})
    
    colorbar
    
    set(gca, 'FontSize',18, ...      % tick labels larger
         'LineWidth',1.5, ...    % axis lines thicker
         'TickDir','out', ...    % ticks outward
         'Box','off')            % remove top/right frame

end
%}

%% kappa and widths influence...

if options.colorGridThetaFixed

    if ~ isscalar(thetas3)
        error('Choose a fixed theta')
    else
        disp(['thetas are ', num2str(theta2), ' ',num2str(thetas3)])
    end

    figure(3)
    %imagesc(D)
    imagesc(widthsAsymRatio,kappas,DFixedTheta)
    ylabel('\kappa', 'Rotation',0)
    
    yticks(kappas)
    
    xlabel('W.A.R.')
    
    xticks(widthsAsymRatio)
    
    %xticklabels({'\pi/12','3\pi/12','5\pi/12'})
    
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
%}
%---------------------
function [h1,h2,h3,th_xi,xi] = loading(kappa, angles, widths, Lx)

[waveName, graphName] = standardNaming(Lx, widths, angles, kappa);

load(graphName,'w','th_xi','th_zeta');
load(waveName,'h')

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



end