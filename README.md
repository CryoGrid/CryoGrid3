# CryoGrid3
CryoGrid 3 is a simple land-surface scheme dedicated to modeling of ground temperatures in permafrost environments.

The model version of is the refernce for the research article "Modelling the degradation of ice-wedges in polygonal tundra under different hydrological conditions" which is intended to be published in the scientific journal "The Cryosphere".

The parameters are set to the default values used in the article. Parameters different from the default values can be specified in the main scirpt "CryoGrid3_xice_mpi" (general parameters) and in the function "get_parallel_variables" (tile-specific parameters) which is part of the module "cryoGridLateral".

To start the program, just run the script "CryoGrid3_xice_mpi". The default output directory is "./runs/".