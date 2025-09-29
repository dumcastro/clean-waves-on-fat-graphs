Fat graph linear wave solver.

This program generates solutions for linear hyperbolic pulses propagating in 2D graph-like shapes using the SC transform.

Domain arguments:
- 3 reach widths
- 2 branching angles
- 1 channel length (all reaches have same length)

Travelling wave argument:
- effective wavelength lambda

The parameter kappa = width-of-reach-1/lambda ties the geometry of the channel with the wave profile.
Length is chosen so that wave fits comfortably within each reach of the channel.

Start by running: runDemo.m
For parameter sweep use: runFullExperiment.m

Directory structure:
│
├── runDemo.m             % Choose a single geometry, create mesh and solve evolution EDP on canonical domain
├── runFullExperiment.m   % Same as runDemo but for a range of geometrical parameters
├── parameter_station.m   % Chooses secondary parameters like number of mesh points and visualization options
├── private/                % Core and auxiliar modular functions
│   ├── createFatGraph.m      % Create *fat graph* object from given geometrical arguments
│   ├── evolveWave.m	      % 'If block' for L shape vs Y shape
│   ├── evolveWave3.m	      % Main time-stepping + finite difference scheme on Y shape. Requires data from createFatGraph
│   ├── evolveWave2.m	      % Same as evolveWave3 but for L shape (sharp angle).
│   ├── processWaveData.m      % Visualization of aspects of wave solutions, plotting
│   ├── processGraphData.m     % Visualization of aspects of graph data, e.g., Jacobian of SC transform
│   ├── outermollif.m     % Creates an epsilon gap in poly to avoid singular points
├── comparison_plot_asym.m   % Compares transmitted waves for asymmetric case

├── GraphData/                % Stores graph data
├── WaveData/                 % Stores wave data
├── External/                 % SC toolbox
├── Export/                   % For videos, figures, plots generated from the data
├── Archive/                  % Some outdated scripts
├── FatGraph.m	          % A class for fat graph objects, with specific lenghts, angles and widths	
├── README.md             % Docs + tutorial


Dependencies:

- Requires Toby Driscoll’s SC Toolbox, included in the External Folder
