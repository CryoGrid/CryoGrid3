% %clear all
% if ~paraFromFile
%     clearvars -except paraFromFile     %added by JAN to enable running with or without config file from same script
% end
close all
profile off
%dbclear if error
%dbstop if error

%import CryoGrid modules (matlab functions)
addpath('modules/cryoGridTechnical/')
addpath('modules/cryoGridInitialize/')
addpath('modules/cryoGridSEB/')
addpath('modules/cryoGridSoil/')
addpath('modules/cryoGridSnow/')
addpath('modules/CryoGridInfiltrationUnfrozenSoil')
addpath('modules/cryoGridExcessIce/')
addpath('modules/cryoGridExcessIceInfiltration')
%addpath('modules/cryoGridRockFields/')
