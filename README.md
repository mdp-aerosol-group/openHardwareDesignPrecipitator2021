# openHardwareDesignPrecipitator2021
Data repository and scripts associated with "Open-Hardware Design and Characterization of an Electrostatic Aerosol Precipitator", Kasparoglu, S. and Petters, M. D., 2021

## Folders
Scripts used to generate the figures: ```src/```

Generated figures: ```src/Figures/```

Data files: ```src/Data/```

Shared files among scripts: ```src/Scripts/```

Computational fluid dynamics simulations: ```zapperFoam/```

## Brief Description
The ```src/``` folder contains a self contained set of scripts for data processing, modeling, and figure written in the julia language. These scripts should be executed from within the folder. The Project.toml and Manifest.toml files summarize the environment (package versions) that was active when running the scripts. Please refer to the julia language documentation for working with environments. The naming convention is ```f0X.jl``` where ```X``` is the figure number.

The ```src/Data/``` folder includes the experimental data from the characterization experiments, stored results from the transfer function model, and fluid velocity data exported from the OpenFoam simulation.

The ```zapperFOAM``` folder contains a case folder with the initialization files for the CFD simulation. Please refer to the OpenFOAM documentation for details. 

To create the mesh call

```zapperFOAM>> blockMesh```

To run the simulation call

```zapperFOAM>> icoFoam```

Velocity data in the folder ```src/Data/``` were exported with ```paraFoam``` for ```t = 4s```.  

Simulations were run using the [docker image](https://hub.docker.com/r/openfoam/openfoam8-paraview56) ```OpenFOAM v8 | Linux``` 




