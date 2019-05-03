% import CryoGrid modules
addpath('modules/cryoGridTechnical/')
addpath('modules/cryoGridInitialize/')
addpath('modules/cryoGridSEB/')
addpath('modules/cryoGridSoil/')
addpath('modules/cryoGridSnow/')
addpath('modules/cryoGridInfiltrationUnfrozenSoil')
addpath('modules/cryoGridExcessIce/')
addpath('modules/cryoGridExcessIceInfiltration')
addpath('modules/cryoGridLateral')
% import Visualization module and load colormap
addpath('modules/cryoGridVisualization')
cm = iLoadColormap( 'cm_blueautumn.mat' );
