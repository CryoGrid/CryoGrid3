clear all;
close all;

%% specify setup file (see subdirectory "./setups/" for further setup files
setupFile = 'SINGLE-TILE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_translational_Nmeso2_xEm0_xHm0_xWm1_Dmeso100_slope0.000_eRes-10.00_setup';

%% run the model
load( setupFile );
CryoGrid3_xice_mpi( SETUP );