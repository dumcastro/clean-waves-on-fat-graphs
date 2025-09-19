function [waveName, graphName] = standardNaming(Lx, widths, angles, kappa)

graphName = ['GraphData/widths=', mat2str(widths),'_angles=', mat2str(round(angles, 3)),'_length=',...
            mat2str(round(Lx, 3)), '.mat'];

waveName = ['WaveData/kappa=', num2str(kappa),'_widths=', mat2str(widths), '_angles=', mat2str(round(angles, 3)), '_length=',...
            mat2str(round(Lx, 3)), '.mat'];

end