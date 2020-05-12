clear all;
close all;

% initialize SETUP struct
SETUP = {};

% meso-topology
SETUP.mesoType = 'translational';                                       % symmetry of the meso-scale tiles
SETUP.Nmeso = 3;                                                        % [-], number of meso-scale tiles
SETUP.Dmeso = 100;                                                      % [m], distance between meso-scale tiles
SETUP.altitude = 20;                                                    % [m], this is the altitude in the central tile of the meso-scale (=center of the island)
SETUP.slope = 0.001;                                                    % [-], slope of the meso-scale gradient

% micro-topology
SETUP.microType = 'polygonCRT'; % 'simple'                              % micro-scale represetnation: 'simple' for homogeneous and 'polygon' for ice-wedge polygons                                 

% climate warming scenario
SETUP.scenario='rcp85'; %'rcp45'

% technical parameters
if strcmp(SETUP.microType , 'polygonCRT')
    SETUP.numRealizations = SETUP.Nmeso .* 3;
else
    SETUP.numRealizations = SETUP.Nmeso;
end
SETUP.syncTimestep=1./24;                                               % [days], length of timestep for lateral exchange
SETUP.startDate = datenum( 1999, 10, 1 );                               % start date of the simulations
SETUP.endDate = datenum( 2099, 12, 31 );                                % end date of the simulations

% switches for lateral transport processes
SETUP.lateral=1;
SETUP.xH=1;
SETUP.xW=1;
SETUP.xS=1;
SETUP.xE=1;
SETUP.xice=1;

% lateral transport on the meso-scale
SETUP.xWmeso=1;
SETUP.xHmeso=0;
SETUP.xEmeso=0;

% snow and hydrology parameters
SETUP.snowDens = 250;                                                   % [kg/m^3], density of snow at deposition
SETUP.relMaxSnow = 0.5;                                                 % [m], maximum height of snowcover relative to the tile with the heighest elevation in the same meso unit
SETUP.fieldCapacity = 0.5;                                              % [-], field capacity of uppermost soil layers

% parameters related to hydrological boundary conditions
SETUP.boundaryConditionFirstMeso ='DarcyReservoirNoInflow';             % this is applied to the last micro unit of the first meso unit
SETUP.eReservoirFirstMeso = -10.0;                                      % this is relative to the initial altitude of the first micro unit of the first meso unit
SETUP.boundaryConditionLastMeso = 'NoBC';                               % this is applied to the last micro unit of the last meso unit
SETUP.eReservoirLastMeso = 0.0;                                         % this is relative to the initial altitude of the first micro unit of the last meso unit

% surface water scheme (meso scale)
SETUP.mesoSurfaceWaterScheme = 'Darcy';                                 % method to calculate surface fluxes at the meso-scale
SETUP.Ksubs=1e-5;                                                       % subsurface hydraulic conductivity in [m/s]
SETUP.Ksurf_meso=1e-2;                                                  % surface hydraulic conductivity (meso-scale) in [m/s]
SETUP.Ksurf_micro=1e-5;                                                 % surface hydraulic conductivity (micro-scale) in [m/s]


if strcmp( SETUP.microType, 'polygonCRT' )
    SETUP.microElev = [ 0.0, 0.2, 0.0 ];                                % [m], relative initial elevations of the center, rim and trough tiles within each meso unit
    SETUP.microArea = 140.*[ 1/3, 1/2, 1/6 ];                           % [m^2], areas of center, rim and trough tiles representing one ice-wedge polygon
else
    SETUP.microElev = [];          
    SETUP.microArea = [];
end

% excess ice placement
SETUP.phiNat = 0.55;
if strcmp(SETUP.microType, 'simple')
    SETUP.Dxice = [ 0.9, 0 ,0  ];                                       % [m], depth of excess ice layer (only first column used)             
    SETUP.thetaIce = [ 0.75, 0, 0 ];                                    % [-], ice content of excess ice layer (only first column used)
elseif strcmp(SETUP.microType, 'polygonCRT')
    SETUP.Dxice = [ 1.0, 0.9, 0.7 ];                                    % [m], depth of excess ice layer (center, rim, trough), relative to the soil surface
    SETUP.thetaIce = [ 0.65, 0.75, 0.95 ];                              % [-], ice content of excess ice layer (center, rim, trough)
end       

% stratigraphy layers [ water/ice content, mineral content, organic content, soily type (1:sand, 2:silt), natural porosity ]
VL     = [ 0.85, 0.00, 0.15, 1, 0.85 ] ;        % vegetation layer
OLrev  = [ 0.75, 0.10, 0.15, 1, 0.75 ] ;        % active layer: organic horizon
MLrev  = [ 0.65, 0.25, 0.10, 2, 0.65 ] ;        % active layer: mineral horizon
ILrev  = [ 0.65, 0.20, 0.15, 1, 0.55 ] ;        % intermediate layer: DTLB, Holocene
TL2    = [ 0.45, 0.50, 0.05, 1, 0.45 ] ;        % older taberite layers
BLrev  = [ 0.10, 0.90, 0.00, 1, 0.10 ] ;        % bedrock

stratigraphyMap = containers.Map( {'CENTER', 'RIM', 'TROUGH', 'SIMPLE' },...
    { [ 0.00, VL; ...
        0.10, OLrev;...
        0.20, MLrev;...
        SETUP.Dxice(1)-0.2, ILrev;...
        SETUP.Dxice(1), [ SETUP.thetaIce(1), round( min( [(1-SETUP.thetaIce(1)) , +0.025+(1-SETUP.thetaIce(1))./2] ), 2 ), round( max( [0, -0.025+(1-SETUP.thetaIce(1))./2] ), 2 ) , 1, SETUP.phiNat ];...
        10.00, TL2;...
        30.00, BLrev ],...
      [ 0.00, VL;...
        0.10, OLrev;...
        0.20, MLrev;...
        SETUP.Dxice(2)-0.2, ILrev;...
        SETUP.Dxice(2), [ SETUP.thetaIce(2), round( min( [(1-SETUP.thetaIce(2)) , +0.025+(1-SETUP.thetaIce(2))./2] ), 2 ), round( max( [0, -0.025+(1-SETUP.thetaIce(2))./2] ), 2 ) , 1, SETUP.phiNat ];...
        10.00, TL2;...
        30.00, BLrev ],...
      [ 0.00, VL; ...
        0.10, OLrev;...
        0.20, MLrev;...
        SETUP.Dxice(3)-0.2, ILrev;...
        SETUP.Dxice(3), [ SETUP.thetaIce(3), round( min( [(1-SETUP.thetaIce(3)) , +0.025+(1-SETUP.thetaIce(3))./2] ), 2 ), round( max( [0, -0.025+(1-SETUP.thetaIce(3))./2] ), 2 ) , 1, SETUP.phiNat ];...
       10.00, TL2;...
       30.00, BLrev ],...
      [ 0.00, VL;...
        0.10, OLrev;...
        0.20, MLrev;...
        SETUP.Dxice(1)-0.2, ILrev;...
        SETUP.Dxice(1), [ SETUP.thetaIce(1), round( min( [(1-SETUP.thetaIce(1)) , +0.025+(1-SETUP.thetaIce(1))./2] ), 2 ), round( max( [0, -0.025+(1-SETUP.thetaIce(1))./2] ), 2 ) , 1, SETUP.phiNat ];...
        10.00, TL2;...
       30.00, BLrev ] ...
    } );

if strcmp(SETUP.microType, 'polygonCRT')
    SETUP.stratigraphy = {  stratigraphyMap([ 'CENTER' ]), ...
                            stratigraphyMap([ 'RIM' ]), ...
                            stratigraphyMap([ 'TROUGH' ]) };
elseif strcmp(SETUP.microType, 'simple')
    SETUP.stratigraphy = {  stratigraphyMap(['SIMPLE']) };
end

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
SETUP.forcingFile = ['Samoylov_' SETUP.scenario '_1901_2300_CryoGrid_windModified_repeat2080to2100until2200.mat'];
[~, SETUP.git_commit_hash] = system('git rev-parse HEAD');

% output directory
SETUP.saveDir = './runs';

% compose runname
SETUP.runName = sprintf( [ 'POLYGON-LANDSCAPE_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' )  '_' SETUP.scenario '_xice%d_xE%d_xH%d_xW%d_xS%d_%s_%s_Nmeso%d_eRes%0.2f' ], ...
        SETUP.xice, SETUP.xE, SETUP.xH, SETUP.xW, SETUP.xS, SETUP.microType, SETUP.mesoType, SETUP.Nmeso, SETUP.eReservoirFirstMeso );

% create ouput directory
mkdir( [ SETUP.saveDir ] );
mkdir( [ SETUP.saveDir '/' SETUP.runName ] );
mkdir( [ SETUP.saveDir '/' SETUP.runName '/settings'] );
mkdir( [ SETUP.saveDir '/' SETUP.runName '/output'] );
mkdir( [ SETUP.saveDir '/' SETUP.runName '/plots'] );
mkdir( [ SETUP.saveDir '/' SETUP.runName '/diagnostics'] );
mkdir( [ SETUP.saveDir '/' SETUP.runName '/finalState'] );

% save setup file
save( [ './setups/' SETUP.runName '_setup.mat' ], 'SETUP'  );
%save( [ SETUP.saveDir '/' SETUP.runName '/' SETUP.runName '_setup.mat' ], 'SETUP'  );