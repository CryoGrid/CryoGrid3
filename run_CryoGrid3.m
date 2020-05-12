clear all;
close all;

%% specify setup file
% homogeneous micro-scale representation
setupFile = './setups/POLYGON-LANDSCAPE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_simple_translational_Nmeso2_eRes-0.10_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"
%setupFile = './setups/POLYGON-LANDSCAPE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_simple_translational_Nmeso2_eRes-10.00_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"
%setupFile = './setups/POLYGON-LANDSCAPE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_simple_translational_Nmeso3_eRes-10.00_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"
% polygon micro-scale representation
%setupFile = './setups/POLYGON-LANDSCAPE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_polygonCRT_translational_Nmeso2_eRes0.00_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"
%setupFile = './setups/POLYGON-LANDSCAPE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_polygonCRT_translational_Nmeso2_eRes-10.00_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"
%setupFile = './setups/POLYGON-LANDSCAPE_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_polygonCRT_translational_Nmeso3_eRes-10.00_setup.mat';        % load temporary setup file created by "setup_CryoGrid3.m"

%% run the model
load( setupFile );
CryoGrid3_xice_mpi( SETUP );