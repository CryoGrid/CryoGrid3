% setup_CryoGrid3.m : This script can be used to create custom setup files for CryoGrid 3 (v1.1.0).
% After running this script, the path and name of the generated setup filed needs to be specified in the script 'run_CryoGrid3.m', which can then be used to start a simulation. 

% clear and close everything before starting
clear all;
close all;

% initialize SETUP struct
SETUP = {};

% technical parameters
SETUP.numRealizations = 3;														% [-], number of coupled tiles (do not change)
SETUP.syncTimestep=1./24;														% [days], timestep for the exchange between laterally coupled tiles
SETUP.startDate = datenum( 1999, 10, 1 );										% start date for the simulation (make sure to initialize correctly)
SETUP.endDate = datenum( 2099, 12, 31);											% end date for the simulation
SETUP.lateral=1;																% switch for exchange of information between laterally coupled tiles (do not change)
SETUP.xH=1;																		% switch for lateral heat transport; set to 0 to disable lateral heat transport
SETUP.xW=1;																		% switch for lateral water transport; set to 0 to disable lateral water transport
SETUP.xS=1;																		% switch for lateral snow redistribution; set to 0 to disable snow redistribution
SETUP.xE=1;																		% switch for lateral sediment transport; set to 0 to diable lateral sediment transport
SETUP.xice=1;																	% switch for the excess ice scheme; set to 0 to ignore presence of excess ground ice

% climate warming scenario
SETUP.scenario='rcp85';															% warming scenario according to representative concentration pathways, possible values: 'rcp26', 'rcp45', 'rcp85' (if changing to 'rcp26' the name of the forcing file (SETUP.forcingFile) needs to be adapted below)

% parameters related to hydrological boundary conditions and snow
SETUP.fieldCapacity = 0.50;														% [-], field capacity, i.e. capacity of the upper soil layers to hold infiltrating water; [Nitzbon et al., 2019, TC]
SETUP.relMaxSnow = 0.40;														% [m], maximum height of snow cover relative to the tile with the heighest elevation
SETUP.relMaxWater = 10.00;														% [m], maximum height of water column above surface that can be simulated (do not change)
SETUP.snowDens = 250;															% [kg/m^3], density of snow at deposition
SETUP.boundaryCondition_T = 'DarcyReservoirNoInflow';							% hydrological boundary condition (do not change)
SETUP.e_Reservoir = -10.0;														% [m], elevation of the external water reservoir relative to the centre tile; 0.0 for water-logged LB or HD, -0.2 for water-logged YD, -10.0 for well-drained LB, HD, YD

% areal fractions, need to add up to 1
SETUP.f_C = 0.3;																% [-], areal fraction of the polygons centres; [Muster et al. (2012), Nitzbon et al., 2019] 
SETUP.f_R = 0.6;																% [-], areal fraction of the polygons rims; [Muster et al. (2012), Nitzbon et al., 2019] 
SETUP.f_T = 0.1;																% [-], areal fraction of the inter-polygon troughs; [Nitzbon et al., 2019] 

% topography
SETUP.e_R = 0.4;																% [m], relative elevation of the rim tile above the centre tile; 0.4 for LB and HD, 0.0 for YD; [Nitzbon et al., 2019, TC]
SETUP.e_T = 0.3;																% [m], relative elevation of the trough tile above the centre tile; 0.3 for LB and HD, 0.0 for YD; [Nitzbon et al., 2019, TC]

% hydraulic conductivities
SETUP.K=1e-5;																	% [m/s], hydraulic conductivity for lateral water fluxes, [Nitzbon et al., 2019, TC]
SETUP.K_subs=SETUP.K;															% [m/s], hydraulic conductivity for lateral subsurface water fluxes, here: no distinction between surface and subsurface
SETUP.K_surf=SETUP.K;															% [m/s], hydraulic conductivity for lateral surface water fluxes, here: no distinction between surface and subsurface
SETUP.K_Reservoir = 2*pi*SETUP.K_subs;											% [m/s], effective hydraulic conductivity to the external water reservoir, [Nitzbon et al., 2019, TC]

% hillslope diffusitivities
SETUP.weight_diffusion = 0;														% [-], weighting factor for diffusive lateral sediment fluxes
SETUP.weight_advection = 1;														% [-], weighting factor for advective lateral sediment fluxes
SETUP.hillslope_diffusivity_land =  3e-10; 										% [m^2/sec]  (approx. 0.01 m^2/yr), subaerial sediment transport coefficient; [ Kessler et al. 2012, JGR ]
SETUP.hillslope_diffusivity_water = 3e-8; 										% [m^2/sec]  (approx. 1.00 m^2/yr), subaqaueous sediment transport coefficient; [ Kessler et al. 2012, JGR ]
SETUP.critical_hillslope_angle = pi/4;											% [rad] criticle angle for which advective fluxes diverge; [ Kessler et al. 2012, JGR ]

% definition of ground stratigraphy layers
VL     = [ 0.85, 0.00, 0.15, 1, 0.85 ] ;        								% vegetation layer
OLrev  = [ 0.75, 0.10, 0.15, 1, 0.75 ] ;        								% active layer: organic horizon
MLrev  = [ 0.65, 0.25, 0.10, 2, 0.65 ] ;        								% active layer: mineral horizon
ILrev  = [ 0.65, 0.20, 0.15, 1, 0.55 ] ;        								% intermediate layer: DTLB, Holocene
ILrevY = [ 0.65, 0.25, 0.10, 1, 0.55 ] ;        								% intermediate layer: Yedoma
IWrev  = [ 0.95, 0.05, 0.00, 1, 0.55 ] ;        								% ice wedge
IWIL21 = [ 0.85, 0.10, 0.05, 1, 0.55 ] ;        								% ice wedge / intermediate mix 2:1
IWIL12 = [ 0.75, 0.15, 0.10, 1, 0.55 ] ;        								% ice wedge / intermediate mix 1:2
TL1    = [ 0.55, 0.35, 0.10, 1, 0.55 ] ;        								% taberite layer
TL2    = [ 0.45, 0.50, 0.05, 1, 0.45 ] ;        								% older taberite layers
BLrev  = [ 0.10, 0.90, 0.00, 1, 0.10 ] ;        								% bedrock

% specification of the ground stratigraphy for the three tiles (CENTER, RIM, TROUGH)
% the first number is the depth in [m] relative to the surface of the tile in which the respective layer starts to extend downwards
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

% initial temperature profile used to initialize the spin-up runs
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
SETUP.forcingFile = ['Samoylov_' SETUP.scenario '_1901_2300_CryoGrid_github.mat'];

% output directory
SETUP.saveDir = './runs';

% compose runname
SETUP.runName = sprintf( [ 'SCENARIO_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' )  '_' SETUP.scenario '_xice%d_xE%d_xH%d_xW%d_xS%d_%s_eRes%0.2f_snowDens%d' ], ...
        SETUP.xice, SETUP.xE, SETUP.xH, SETUP.xW, SETUP.xS, SETUP.boundaryCondition_T, SETUP.e_Reservoir, SETUP.snowDens );

% save setup file
save( [ './setups/' SETUP.runName '_setup.mat' ], 'SETUP'  );
