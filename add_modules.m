<<<<<<< HEAD
% %clear all
% if ~paraFromFile
%     clearvars -except paraFromFile     % enable running with or without config file from same script
% end

clear all
close all

profile off


%import CryoGrid modules (matlab functions)
=======
clear all
close all

% import CryoGrid modules
>>>>>>> origin/xice_mpi_polygon_TC
addpath('modules/cryoGridTechnical/')
addpath('modules/cryoGridInitialize/')
addpath('modules/cryoGridSEB/')
addpath('modules/cryoGridSoil/')
addpath('modules/cryoGridSnow/')
<<<<<<< HEAD
addpath('modules/CryoGridInfiltrationUnfrozenSoil')
=======
addpath('modules/cryoGridInfiltrationUnfrozenSoil')
>>>>>>> origin/xice_mpi_polygon_TC
addpath('modules/cryoGridExcessIce/')
addpath('modules/cryoGridExcessIceInfiltration')
addpath('modules/cryoGridLateral')
