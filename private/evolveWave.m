function [H,h] = evolveWave(kappa, Lx, widths, angles, options)

if numel(widths) == 3
    [H,h] = evolveWave3(kappa, Lx, widths, angles,options);
elseif numel(widths) == 2
    [H,h] = evolveWave2(kappa, Lx, widths, angles,options);
end

end