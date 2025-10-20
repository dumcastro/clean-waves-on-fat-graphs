
travelDistance = 12;

%% Parameter sweep view options
sweepVisOptions = struct();
sweepVisOptions.colorGridWidthsFixed = false; 
sweepVisOptions.colorGridThetaFixed = false;
sweepVisOptions.colorGridKappaFixed = false;

sweepVisOptions.deltaHeightPlot = false; % (for fixed widths only)

%% Wave parameters
wave_options = struct();
wave_options.T = 200; %Final time of execution
wave_options.want_save = true; 

%% Wave view options
wave_vis_options = struct();
wave_vis_options.play_movie_phys = true;
wave_vis_options.play_movie_canonical = false;
wave_vis_options.jmp_xi = 1;  % (e.g. if this is 2 the program plots every other xi data point)
wave_vis_options.jmp_zeta = 1; % (same as above but for zeta)
wave_vis_options.twoD_plot = false;
wave_vis_options.twoD_animation = false;

%% Graph parameters
graph_options = struct();
graph_options.ep = 0.01;
graph_options.Nzeta = 40; 
% Nzeta: number of mesh points along the transver direction
% this gives dzeta and we set dxi = dzeta

%% Graph view options
graph_vis_options = struct();
graph_vis_options.jmp_xi = wave_vis_options.jmp_xi;
graph_vis_options.jmp_zeta = wave_vis_options.jmp_zeta;
options.plotJ = true;
options.visualizeGrid = true;
