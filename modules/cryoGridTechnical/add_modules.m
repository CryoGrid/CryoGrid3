clear all
close all
profile off
dbclear if error
%dbstop if error

%import CryoGrid modules (matlab functions)
addpath('modules/cryoGridTechnical/')
addpath('modules/cryoGridInitialize/')
addpath('modules/cryoGridSEB/')
addpath('modules/cryoGridSoil/')
addpath('modules/cryoGridSnow/')

addpath('modules/cryoGridExcessIce/')
addpath('modules/cryoGridRockFields/')
