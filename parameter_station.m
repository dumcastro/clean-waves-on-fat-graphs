%% Parameter wall
%% Wave parameters
wave_options = struct();
wave_options.point_source = false;
wave_options.T = 400; %Final time of execution
wave_options.want_save = true;

%% Wave view options
wave_vis_options = struct();
wave_vis_options.play_movie_phys = false;
wave_vis_options.play_movie_canonical = false;
wave_vis_options.jmp_xi = 2;
wave_vis_options.jmp_zeta = 4;
wave_vis_options.twoD_plot = true;
wave_vis_options.twoD_animation = false;


%% Graph parameters
graph_options = struct();
graph_options.ep = 0.2;
%graph_options.Nzeta = widths(1)*4;
graph_options.Nzeta = 40;
%graph_options.dxi = widths(1)*0.04;
%graph_options.dzeta = widths(1)*0.04;
graph_options.plot_flag = false;

%% Graph view options
graph_vis_options = struct();
graph_vis_options.jmp_xi = wave_vis_options.jmp_xi;
graph_vis_options.jmp_zeta = wave_vis_options.jmp_zeta;
options.plotJ = true;
options.visualizeGrid = true;

%% Comparison plot options (asym)
parameterSweep = false;
deltaHeightPlot = true;

%% Extra
travel_distance = 10;
