# CryoGrid3
CryoGrid 3 is a simple land-surface scheme dedicated to modeling of ground temperatures in permafrost environments.

The model version of this branch ("xice_mpi_polygon_scenarios") has been used for the simulations related to a research article on the response of ice-rich permafrost to a warming climate.

The forcing data are contained in the directory "./forcing".

The settings of the runs presented in the research article are provided in the directory "./setups".

The script "setup_CryoGrid3.m" can be used to create further settings files which will be saved to "./setups" by default.

The script "run_CryoGrid3.m" can be used to run the model. For this, the desired settings file needs to be specified in the script. 

The default output directory is "./runs/".