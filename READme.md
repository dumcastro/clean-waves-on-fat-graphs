Fat graph linear wave solver.

This program generates solutions for linear hyperbolic pulses propagating in "fat graphs".

Start by running: runDemo.m

Directory structure:
│
├── README.md                 % Docs + tutorial
├── Scripts/                  % Top-level scripts
│   ├── runDemo.m             % Runs each function on demmand, from geometrical arguments
│   ├── runFullExperiment.m   % Creates a fatgraph object, executes evolution and visualize data, from geometrical arguments
│   ├── parameter_station.m   % Chooses secondary parameters
│
├── Functions/                % Core modular functions
│   ├── FatGraph/	      % A class for fat graph objects, with specific lenghts, angles and widths	
│   ├── evolveWave/	      % Solve linear hyperbolic PDE for a given initial pulse on one of the branches
│   ├── createFatGraph/       % Create fat graph object from given geometrical arguments
│   ├── processWaveData/      % Visualization of aspects of wave solutions
│   ├── processGraphData/     % Visualization of aspects of graph data
│
├── GraphData/                % Stores graph data
├── WaveData/                 % Stores wave data
├── External/                 % SC toolbox
├── Export/                   % For videos, figures, plots generated from the data


Dependencies:

- Requires Toby Driscoll’s SC Toolbox, included in the External Folder
