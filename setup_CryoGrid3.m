clear all;
close all;

% initialize SETUP struct
SETUP = {};

% technical parameters
SETUP.numRealizations = 3;
SETUP.syncTimestep=1./24;
SETUP.startDate = datenum( 1999, 10, 1 );
SETUP.endDate = datenum( 2099, 12, 31);
SETUP.lateral=1;
SETUP.xH=1;
SETUP.xW=1;
SETUP.xS=1;
SETUP.xE=1;
SETUP.xice=1;

% climate warming scenario
SETUP.scenario='rcp85';

% parameters related to hydrological boundary conditions and snow
SETUP.fieldCapacity = 0.50;
SETUP.relMaxSnow = 0.40;
SETUP.relMaxWater = 10.00;
SETUP.snowDens = 250;
SETUP.boundaryCondition_T = 'DarcyReservoirNoInflow';
SETUP.e_Reservoir = -10.0;

% areal fractions
SETUP.f_C = 0.3;
SETUP.f_T = 0.1;
SETUP.f_R = 1.0-SETUP.f_T-SETUP.f_C;

% topography
SETUP.e_R = 0.4;
SETUP.e_T = 0.3;

% hydraulic conductivities
SETUP.K=1e-5;
SETUP.K_subs=1e-5;
SETUP.K_surf=1e-5;
SETUP.K_Reservoir = 2*pi*SETUP.K_subs;

% hillslope diffusitivities
SETUP.weight_diffusion = 0;
SETUP.weight_advection = 1;
SETUP.hillslope_diffusivity_land =  3e-10; % [m^2/sec] 3e-10 m^2/sec approx. 0.01 m^2/yr, reference: [ Kessler et al. 2012, JGR ]
SETUP.hillslope_diffusivity_water = 3e-8; % [m^2/sec]  3e-8  m^2/sec approx  1.00 m^2/yr, reference: [ Kessler et al. 2012, JGR ]
SETUP.critical_hillslope_angle = pi/4;

% stratigraphy layers
VL     = [ 0.85, 0.00, 0.15, 1, 0.85 ] ;        % vegetation layer
OLrev  = [ 0.75, 0.10, 0.15, 1, 0.75 ] ;        % active layer: organic horizon
MLrev  = [ 0.65, 0.25, 0.10, 2, 0.65 ] ;        % active layer: mineral horizon
MLrevSand  = [ 0.65, 0.25, 0.10, 1, 0.65 ] ;    % active layer: mineral horizon
ILrev  = [ 0.65, 0.20, 0.15, 1, 0.55 ] ;        % intermediate layer: DTLB, Holocene
ILrevY = [ 0.65, 0.25, 0.10, 1, 0.55 ] ;        % intermediate layer: Yedoma
IWrev  = [ 0.95, 0.05, 0.00, 1, 0.55 ] ;        % ice wedge
IWIL21 = [ 0.85, 0.10, 0.05, 1, 0.55 ] ;        % ice wedge / intermediate mix 2:1
IWIL12 = [ 0.75, 0.15, 0.10, 1, 0.55 ] ;        % ice wedge / intermediate mix 1:2
TL1    = [ 0.55, 0.35, 0.10, 1, 0.55 ] ;        % taberite layer
TL2    = [ 0.45, 0.50, 0.05, 1, 0.45 ] ;        % older taberite layers
BLrev  = [ 0.10, 0.90, 0.00, 1, 0.10 ] ;        % bedrock

stratigraphyMap = containers.Map( {'CENTER', 'RIM', 'TROUGH' },...
    { [ 0.00, VL; ...
        0.10, OLrev;...
        0.20, MLrev;...
        0.60, ILrevY;...
       20.00, TL1;...
       30.00, TL2;...
       50.00, BLrev ],...
      [ 0.00, VL; ...
        0.10, OLrev;...
        0.20, MLrev;...
        0.60, ILrevY;...
        1.00, IWIL21;...
       20.00, TL1;...
       30.00, TL2;...
       50.00, BLrev ],...
      [ 0.00, VL; ...
        0.10, OLrev;...
        0.20, MLrev;...
        0.60, ILrevY;...
        1.00, IWrev;...         
       20.00, TL1;...
       30.00, TL2;...
       50.00, BLrev ],...
    } );

SETUP.stratigraphy = { stratigraphyMap([ 'CENTER' ]), ...
    stratigraphyMap([ 'RIM' ]), ...
    stratigraphyMap([ 'TROUGH' ]) };

SETUP.Tinitial = {};
for i=1:SETUP.numRealizations
    SETUP.Tinitial{i} = [   -2     5    ;...
                             0     0    ;...
                             2    -2    ;...
                             5    -7    ;...
                            10    -9    ;...
                            25    -9    ;...
                           100    -8    ;...
                          1100    10.2  ];
end

% forcing data 
SETUP.forcingFile = ['Samoylov_' SETUP.scenario '_1901_2300_CryoGrid_windModified_repeat2090sTo2160.mat'];

% output directory
SETUP.saveDir = './runs';

% compose runname
SETUP.runName = sprintf( [ 'SCENARIO_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' )  '_' SETUP.scenario '_xice%d_xE%d_xH%d_xW%d_xS%d_%s_eRes%0.2f_snowDens%d' ], ...
        SETUP.xice, SETUP.xE, SETUP.xH, SETUP.xW, SETUP.xS, SETUP.boundaryCondition_T, SETUP.e_Reservoir, SETUP.snowDens );

% save setup file
save( [ './setups/' SETUP.runName '_setup.mat' ], 'SETUP'  );