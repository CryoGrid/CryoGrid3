function PARA = get_parallel_variables(PARA, SETUP)

index = labindex;

% auxiliary calculations for geometry
% input: typical total polygon area + areal weights (as integers) of center/rim/troughs
area_tot = 140.0; % typical area of polygon center [Cresto Aleina / Muster ]
PARA.ensemble.weight = [SETUP.f_C, SETUP.f_R, SETUP.f_T].*100;

% areas
PARA.ensemble.area = PARA.ensemble.weight ./ sum(PARA.ensemble.weight) .* area_tot ; % in m^2
area_C = PARA.ensemble.area(1);
area_R = PARA.ensemble.area(2);
area_T = PARA.ensemble.area(3);

% distances
distance_CR = (0.5.*area_C + 0.25.*area_R) ./ sqrt( area_tot ); % in m
distance_RT = (0.5.*area_T + 0.25.*area_R) ./ sqrt( area_tot );
% hydraulic distances
halfWidth_R = (0.25.*area_R) ./ sqrt( area_tot );

% perimeters % assuming hexagonal shapes of centers and polygons
perimeter_CR = 6 .* sqrt( 2 .* area_C ./ (3 .* sqrt(3) ) );   % assuming hexagonal shape of center/rim interface %2 .* pi .* ( diameter_C ./2);    % assuming circular shape of polygon centers
perimeter_RT = 6. * sqrt( 2 .* (area_C+area_R) ./ (3 .* sqrt(3) ) ); % assuming hexagonal shape of rim/trough interface

% geometric relations
PARA.ensemble.distanceBetweenPoints = [ 0, distance_CR, 0; distance_CR, 0, distance_RT; 0, distance_RT, 0 ]; %diameter .* ( ones(numlabs) - eye(numlabs) );% [0, 4, 4; 4, 0, 4 ; 4, 4, 0];  %in m; put 0 for all non-connected ensemble members
A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)

% topographical relations
altitude_C = 20.0;
elevation_R = 0.4;
elevation_T = 0.3;
altitude_R = altitude_C + elevation_R;
altitude_T = altitude_C + elevation_T;

PARA.ensemble.initial_altitude = [ altitude_C, altitude_R, altitude_T ]; %[20.0, 21.0, 20.5]; %in m a.s.l., this is the reference for the "zero" position of the grids
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to heat exchange
PARA.ensemble.thermal_contact_length = [0, perimeter_CR, 0; perimeter_CR, 0, perimeter_RT; 0, perimeter_RT, 0 ]; % [ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;

% parameters related to water exchange
PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity .* A;%[ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %in m/sec % [Roth: 1e-5 for saturated silt, 		2.2e-5 for saturated sand]
PARA.ensemble.water_table_altitude = [NaN NaN NaN];  %initialize somehow;

%PARA.ensemble.max_water_flux= [0 0];   %in m water equivalent
PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;%[ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ];
PARA.ensemble.active_layer_depth_altitude = [NaN NaN NaN];
%PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;
PARA.ensemble.hydraulicDistance = [ 0, halfWidth_R, 0; halfWidth_R, 0, halfWidth_R; 0, halfWidth_R, 0 ];

% parameters related to snow exchange
%PARA.ensemble.snow_diffusivity = PARA.snow.diffusivity;
%PARA.ensemble.relative_max_snow_height = 0.2;
PARA.ensemble.immobile_snow_height = [0.1, 0.1, 0.1 ]; %, 0.2 ];  %in m %this replaces PARA.snow.maxSnow ?
PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);
PARA.ensemble.snow_contact_length = PARA.ensemble.thermal_contact_length;

% parameters related to infiltration scheme
PARA.ensemble.external_water_flux=[0, 0, SETUP.extFlux.*area_tot./area_T ] ;  %in m/day
PARA.ensemble.rootDepth = [ 0.2, 0.1, 0.2 ];
PARA.ensemble.fieldCapacity = [ 0.5, 0.5, 0.5 ];

% location-specific fix parameter values
PARA.location.initial_altitude = PARA.ensemble.initial_altitude(index);
PARA.soil.externalWaterFlux = PARA.ensemble.external_water_flux(index);
PARA.soil.rootDepth = PARA.ensemble.rootDepth(index);
PARA.soil.fieldCapacity= PARA.ensemble.fieldCapacity(index);
% location-specific dynamic auxiliary variables
PARA.location.area = PARA.ensemble.area(index);
PARA.location.altitude = PARA.ensemble.altitude(index);
PARA.location.surface_altitude = PARA.ensemble.surface_altitude(index);
PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(index);
PARA.location.active_layer_depth_altitude = PARA.ensemble.active_layer_depth_altitude(index);
% location-specific dynamic common thresholds
PARA.location.absolute_maxWater_altitude = [max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater];
PARA.location.absolute_maxSnow_altitude = [max( PARA.ensemble.altitude ) + PARA.snow.relative_maxSnow];


% different stratigraphies
depth_xice_C = SETUP.d_xice_C;
vwc_xice_C = 0.65;
depth_xice_R = SETUP.d_xice_R;
vwc_xice_R = 0.75;
depth_xice_T = SETUP.d_xice_T;
vwc_xice_T = 0.90;
natPor = 0.5;
stratigraphyMap= containers.Map( { 'DEFAULT', 'WET', 'DRY', 'CENTER', 'RIM', 'TROUGH' }, ...
    { [ 0.0   0.40    0.10    0.15    1   0.75;...
    0.15  0.65    0.30    0.05    1   0.65;...
    0.9   0.65    0.30    0.05    1   0.65;...
    9.0   0.30    0.70    0.00    1   0.30     ], ...
    [ 0.0   0.85    0.00    0.15    1   0.85;...
    0.15  0.75    0.2     0.05    1   0.75;...
    0.30  0.65    0.3     0.05    2   0.65;...
    0.9   0.65    0.3     0.05    1   0.65;...
    9.0   0.30    0.70    0.00    1   0.30     ], ...
    [ 0.0   0.50    0.10    0.15    1   0.75;...
    0.1   0.65    0.30    0.05    2   0.65;...
    0.9   0.65    0.30    0.05    1   0.65;...
    9.0   0.30    0.70    0.00    1   0.30     ],...
    [ 0.0             0.85            0.00    0.15    1   0.85;...
    0.15            0.75            0.20    0.05    1   0.75;...
    0.30            0.65            0.30    0.05    2   0.65;...
    depth_xice_C    vwc_xice_C      0.30    0.05    1   natPor;...
    9.0             0.30            0.70    0.00    1   0.30     ], ...
    [ 0.0             0.50            0.10    0.15    1   0.75;...
    0.1             0.65            0.30    0.05    2   0.65;...
    depth_xice_R    vwc_xice_R      0.20    0.05    1   natPor;...
    9.0+elevation_R 0.30            0.70    0.00    1   0.30     ], ...
    [ 0.0             0.40            0.00    0.15    1   0.85;...
    depth_xice_T    vwc_xice_T      0.05    0.05    1   natPor;...
    9.0+elevation_T 0.30            0.70    0.00    1   0.30     ] } );

PARA.soil.layer_properties = { stratigraphyMap('CENTER'), ...
    stratigraphyMap('RIM'), ...
    stratigraphyMap('TROUGH') };

PARA.soil.layer_properties = PARA.soil.layer_properties{index};

% typical profile for beginning of October
PARA.Tinitial = [  -2     5   ;...
    0     0   ;...
    2    -2   ;...
    5    -7   ;...
    10    -9  ;...
    25    -9   ;...
    100    -8   ;...
    1100    10.2   ];      % the geothermal gradient for Qgeo=0.05W/m² and K=2.746W/Km is about 18.2 K/km


%     % different initial conditions
%      PARA.Tinitial = [-5     5    5     5;...
%                        0    -5   -5     -5;...
%                        1    -5   -5     -5;...
%                       10    -8   -8     -8;...
%                       20   -10  -10     -10;...
%                      100   -10  -10     -10;...
%                     2000    10   10     10];
%
%      PARA.Tinitial=[PARA.Tinitial(:,1) PARA.Tinitial(:, 1+index)];

end
