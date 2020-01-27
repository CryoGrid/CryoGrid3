
clear all;
close all;

%% run the model
setupFile = './setups/SCENARIO_REV6b_194910-211912_rcp45_xice0_xE0_xH0_xW0_xS0_reference_snowDens250_maxSnow1.00_maxWater0.00_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"
load( setupFile );
CryoGrid3_xice_mpi( SETUP );
