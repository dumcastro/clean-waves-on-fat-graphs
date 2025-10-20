function [] = parameterSweep(widthTriplets, kappas, theta2, thetas3)

anglesStart = [0, pi - theta2, pi + 0];
angles = anglesStart;

parameterStation

for kk = 1:length(widthTriplets)
    widths = widthTriplets{kk};
    for ii = 1:length(kappas)
        kappa = kappas(ii);
        lambda_f = widths(1)/kappa;
        Lx = lambda_f * (travelDistance + 1) / 2;
    
        for jj = 1:length(thetas3)
            angles(3) = anglesStart(3) + thetas3(jj);
    
            createFatGraph(Lx, widths, angles,graph_options);
    
            evolveWave(kappa, Lx, widths, angles,wave_options);
    
            close all
        end
        angles = anglesStart;

    end
end



end