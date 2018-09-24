% in this function all parameters which are different for the different
% tiles, are specified and overwrite previously set values; in addtion, the
% topological relations (in this case for polygonal tundra) are specified
function PARA = get_parallel_variables(PARA)

index = labindex;

% topological relations (areas, distances, contact lenghts)
% input: typical total polygon area + areal fractions of center/rim/troughs
area_tot = 140.0;   % typical area of polygon center [Cresto Aleina / Muster ]
f_C = 0.3;          % areal fraction of the centres
f_R = 0.6;          % areal fraction of the rims
f_T = 0.1;          % areal fraction of the troughs
PARA.ensemble.weight = round( [f_C, f_R, f_T], 1 ) .* 100;  % make sure to have integers as weights

PARA.ensemble.area = PARA.ensemble.weight ./ sum(PARA.ensemble.weight) .* area_tot ; % in m^2
area_C = PARA.ensemble.area(1);
area_R = PARA.ensemble.area(2);
area_T = PARA.ensemble.area(3);

% thermal distances
distance_CR = (0.5.*area_C + 0.25.*area_R) ./ sqrt( area_tot ); % in m
distance_RT = (0.5.*area_T + 0.25.*area_R) ./ sqrt( area_tot );
% hydraulic distances
halfWidth_R = (0.25.*area_R) ./ sqrt( area_tot );
PARA.ensemble.distanceBetweenPoints = [ 0, distance_CR, 0; distance_CR, 0, distance_RT; 0, distance_RT, 0 ];
A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)

% perimeters (=contact lengths);
perimeter_CR = 6 .* sqrt( 2 .* area_C ./ (3 .* sqrt(3) ) );             % assuming hexagonal shape of center/rim interface
perimeter_RT = 6. * sqrt( 2 .* (area_C+area_R) ./ (3 .* sqrt(3) ) );    % assuming hexagonal shape of rim/trough interface

% micro-topography
altitude_C = 20.0;
altitude_R = 20.4;
altitude_T = 20.3;
elevation_res = 0.0;  % elevation of the external water reservoir, relative to the altitude of the center, i.e. altitude_res=altitude_C+elevation_res

PARA.ensemble.initial_altitude = [ altitude_C, altitude_R, altitude_T ]; %in m a.s.l., this is the reference for the "zero" position of the grids
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to HEAT exchange
PARA.ensemble.thermal_contact_length = [0, perimeter_CR, 0; perimeter_CR, 0, perimeter_RT; 0, perimeter_RT, 0 ];
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;
PARA.ensemble.thermal_contact_DepthOverDistance = 0; % this factor is used to control the depth below which heat exchange takes place: contact_altitude = min(altitudes) - DepthOverDistance * Distance_ij

% parameters related to WATER exchange
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * A;
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
PARA.ensemble.infiltration_altitude = nan(1, numlabs);
PARA.ensemble.hydraulicDistance = [ 0, halfWidth_R, 0; halfWidth_R, 0, halfWidth_R; 0, halfWidth_R, 0 ];

K_res = 5e-5;
boundaryCondition={'NoBC','NoBC', 'DarcyReservoir'};
Darcy_elevation= [ nan nan altitude_C + elevation_res ];        % Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers withour this boundary condition
Darcy_fluxFactor=[ nan nan K_res ];                 % Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance    Defined for now like this, lets see if we wantto define it differently
PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
    end
end

% parameters related to SNOW exchange
PARA.ensemble.immobile_snow_height = [ 0.1, 0.1, 0.1 ]; %in [m]
PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);

% parameters related to infiltration scheme
PARA.ensemble.external_water_flux=[0, 0, 0 ];
PARA.ensemble.rootDepth = [ 0.2, 0.1, 0.2 ];
PARA.ensemble.fieldCapacity = [ PARA.soil.fieldCapacity, PARA.soil.fieldCapacity, PARA.soil.fieldCapacity ];

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
PARA.location.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
PARA.location.soil_altitude = PARA.ensemble.soil_altitude(index);
% location-specific dynamic common thresholds
PARA.location.absolute_maxWater_altitude = [max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater];
PARA.location.absolute_maxSnow_altitude = [max( PARA.ensemble.altitude ) + PARA.snow.relative_maxSnow];

% tile-specific stratigraphies
natPor = 0.55;
stratigraphyMap= containers.Map( { 'CENTER', 'RIM', 'TROUGH' }, ...
    { [ 0.0             0.85            0.00    			0.15    1   0.85;...
        0.15            0.75            0.20    			0.05    1   0.75;...
        0.30            0.65            0.30    			0.05    2   0.65;...
        0.90            0.65            0.30                0.05    1   natPor;...
        9.0             0.30            0.70    			0.00    1   0.30     ], ...
        [ 0.0           0.50            0.10    			0.15    1   0.75;...
        0.10            0.65            0.30    			0.05    2   0.65;...
        0.60            0.75            0.20                0.05    1   natPor;...
        9.40            0.30            0.70    			0.00    1   0.30     ], ...
        [ 0.0           0.50            0.00    			0.15    1   0.85;...
        0.2             0.75            0.20                0.05	1   natPor; ...
        0.5             0.90            0.05                0.05    1   natPor;...
        9.30            0.30            0.70    			0.00    1   0.30     ] } );

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
end
