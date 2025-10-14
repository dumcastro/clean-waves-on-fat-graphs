%% Parameter wall

travelDistance = 10;

%% Wave parameters
wave_options = struct();
wave_options.T = 100; %Final time of execution
wave_options.want_save = true;

%% Wave view options
wave_vis_options = struct();
wave_vis_options.play_movie_phys = true;
wave_vis_options.play_movie_canonical = false;
wave_vis_options.jmp_xi = 4;
wave_vis_options.jmp_zeta = 4;
wave_vis_options.twoD_plot = false;
wave_vis_options.twoD_animation = false;

%% Graph parameters
graph_options = struct();
graph_options.ep = 0.01;
graph_options.Nzeta = 30;

%% Graph view options
graph_vis_options = struct();
graph_vis_options.jmp_xi = wave_vis_options.jmp_xi;
graph_vis_options.jmp_zeta = wave_vis_options.jmp_zeta;
options.plotJ = true;
options.visualizeGrid = true;

%% Comparison plot options (asym)
parameterSweep = true;
deltaHeightPlot = false;
