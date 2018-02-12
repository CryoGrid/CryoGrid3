% %clear all
% if ~paraFromFile
%     clearvars -except paraFromFile     %added by JAN to enable running with or without config file from same script
% end

clear all
close all

profile off


%import CryoGrid modules (matlab functions)
addpath('modules/cryoGridTechnical/')
addpath('modules/cryoGridInitialize/')
addpath('modules/cryoGridSEB/')
addpath('modules/cryoGridSoil/')
addpath('modules/cryoGridSnow/')
addpath('modules/CryoGridInfiltrationUnfrozenSoil/')
%tsvd new modules added
addpath('modules/cryoGridExcessIce/')
addpath('modules/cryoGridExcessIceInfiltration')
addpath('modules/cryoGridLateral')
%tsvd
addpath('modules/cryoGridLake/')  