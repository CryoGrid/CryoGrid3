# CryoGrid3
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3628658.svg)](https://doi.org/10.5281/zenodo.3628658)

## Background

**CryoGrid 3** is a land-surface scheme dedicated to modeling of ground temperatures in permafrost environments. Its excess ice module **Xice** is capable of simulating ground subsidence and thermokarst lake formation due to melting of excess ground ice. The basic version of **CryoGrid 3** and **Xice** is described in the following article:

Westermann, S., Langer, M., Boike, J., Heikenfeld, M., Peter, M., Etzelmüller, B., & Krinner, G. (2016). Simulating the thermal regime and thaw processes of ice-rich permafrost ground with the land-surface model CryoGrid 3. *Geosci. Model Dev.*, 9(2), 523–546. [https://doi.org/10.5194/gmd-9-523-2016](https://doi.org/10.5194/gmd-9-523-2016)

Version `v1.0.0` has been extended by a hydrology scheme for unfrozen ground conditions, and schemes for the lateral transport of heat, water, and snow between adjacent parts of the simulated environment. It has been set-up to study the degradation of ice-wedge polygons as described in the following article:

Nitzbon, J., Langer, M., Westermann, S., Martin, L., Aas, K. S., & Boike, J. (2019). Pathways of ice-wedge degradation in polygonal tundra under different hydrological conditions. *The Cryosphere*, 13, 1089–1123. [https://doi.org/10.5194/tc-13-1089-2019](https://doi.org/10.5194/tc-13-1089-2019)

The version of this release `v1.1.0` has been extended by a scheme for the lateral transport of sediment within the simulated environment. It has been used for the simulations related to a research article on the response of ice-rich permafrost to a warming climate.


## Instructions

### System requirements
- The code needs to be executed using the software *MATLAB* with the *Statistics and Machine Learning Toolbox* and the *Parallel Computing Toolbox*.
- The code has been tested using the *MATLAB* versions `R2017a` and `R2018b`.
- To run the code, a CPU with at least 4 cores is required. 

### Installation
A dedicated installation or compilation of the code is not necessary.

### Usage
- The forcing data needed to run the model are contained in the directory `./forcing`.
- The script `setup_CryoGrid3.m` can be used to create setup files which will be saved to the directory `./setups` by default. In that script different model parameters can be varied.
- The script `run_CryoGrid3.m` can be used to run the model code. For this, the desired setup file needs to be specified in the script. Then the script can be executed.
- The default output directory is `./runs`.

### Reproduction of results
- The setup files for the runs presented in the research article are provided in the directory `./setups`.
- To reproduce the results, the respective setup files need to be specified in the script `run_CryoGrid3.m` which then needs to be executed.
- Note that the usual computation time for a 100-year run is about 2 to 3 days on a state-of-the-art CPU.
