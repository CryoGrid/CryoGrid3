% in this function all parameters which are different for the different
% tiles, are specified and overwrite previously set values; in addtion, the
% topological relations are specified
function PARA = get_parallel_variables(PARA, SETUP)

index = labindex;

% topological relations
area_tot = 140.0;
f_C = SETUP.f_C;
f_R = SETUP.f_R;
f_T = SETUP.f_T;
PARA.ensemble.weight = round( [f_C, f_R, f_T], 1 ) .* 100;  % make sure to have integers as weights
PARA.ensemble.area = PARA.ensemble.weight ./ sum(PARA.ensemble.weight) .* area_tot ; % in m^2
area_C = PARA.ensemble.area(1);
area_R = PARA.ensemble.area(2);
area_T = PARA.ensemble.area(3);

% polygon geometry
if SETUP.polygon_geometry == 1 % % HEXAGONAL GEOMETRY

    % thermal distances
    distance_CR = (0.5.*area_C + 0.25.*area_R) ./ sqrt( area_tot ); % in m
    distance_RT = (0.5.*area_T + 0.25.*area_R) ./ sqrt( area_tot );
    % hydraulic distances
    halfWidth_R = (0.25.*area_R) ./ sqrt( area_tot );
    PARA.ensemble.distanceBetweenPoints = [ 0, distance_CR, 0; distance_CR, 0, distance_RT; 0, distance_RT, 0 ];

    % perimeters % assuming hexagonal shapes of centers and polygons
    perimeter_CR = 6 .* sqrt( 2 .* area_C ./ (3 .* sqrt(3) ) );   % assuming hexagonal shape of center/rim interface %2 .* pi .* ( diameter_C ./2);    % assuming circular shape of polygon centers
    perimeter_RT = 6. * sqrt( 2 .* (area_C+area_R) ./ (3 .* sqrt(3) ) ); % assuming hexagonal shape of rim/trough interface

    % parameters related to HEAT exchange
    PARA.ensemble.thermal_contact_length = [0, perimeter_CR, 0; perimeter_CR, 0, perimeter_RT; 0, perimeter_RT, 0 ]; % [ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %
    PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;

    % parameters related to WATER exchange
    PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
    PARA.ensemble.hydraulicDistance = [ 0, halfWidth_R, 0; halfWidth_R, 0, halfWidth_R; 0, halfWidth_R, 0 ];

elseif SETUP.polygon_geometry == 2 % CIRCULAR CROSS-SECTION

    radius_C = sqrt( area_C / pi );
    radius_R = sqrt( (area_C+area_R) / pi );
    radius_T = sqrt( (area_C+area_R+area_T) / pi );

    perimeter_CR = 2*pi*radius_C;
    perimeter_RT = 2*pi*radius_R;

    distance_CR = ( radius_C + radius_R ) / 2;
    distance_RT = radius_T - radius_R/2 - radius_C/2;

    % distances
    PARA.ensemble.distanceBetweenPoints = [ 0, distance_CR, 0; distance_CR, 0, distance_RT; 0, distance_RT, 0 ];
    PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;
    PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;

    % contact lengths
    PARA.ensemble.thermal_contact_length = [0, perimeter_CR, 0; perimeter_CR, 0, perimeter_RT; 0, perimeter_RT, 0 ]; % [ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %
    PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
end

% topographical relations
A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)

altitude_C = 0.0; %20.0;
elevation_R = SETUP.e_R;
elevation_T = SETUP.e_T;
altitude_R = altitude_C + elevation_R;
altitude_T = altitude_C + elevation_T;
elevation_Reservoir = SETUP.e_Reservoir;

PARA.ensemble.initial_altitude = [ altitude_C, altitude_R, altitude_T ]; %in m a.s.l., this is the reference for the "zero" position of the grids
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to WATER exchange
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * A;
PARA.ensemble.hydraulic_conductivity_subs = PARA.soil.hydraulic_conductivity_subs * A;
PARA.ensemble.hydraulic_conductivity_surf = PARA.soil.hydraulic_conductivity_surf * A;
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.infiltration_altitude = nan(1, numlabs);

boundaryCondition={'NoBC','NoBC', SETUP.boundaryCondition_T}; 		% set to 'DarcyReservoir' for an external water reservoir
Darcy_elevation=[ nan, nan, altitude_C+elevation_Reservoir ] ; 				% Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers without this boundary condition
Darcy_fluxFactor= [ nan nan SETUP.K_Reservoir ];		% Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance. Defined for now like this, lets see if we wantto define it differently. NaN for workers without this boundary condition

%boundaryCondition={'NoBC','NoBC', 'DarcyReservoirDynamicConductivity'}; 		% set to 'DarcyReservoir' for an external water reservoir
% Darcy_fluxFactorMax= [ nan nan 2*pi*PARA.soil.hydraulic_conductivity_surf ];    % factor 2pi from assuming ciruclar reservoir in distance D
% Darcy_fluxFactorMin=[ nan nan 2*pi*PARA.soil.hydraulic_conductivity_subs ];
% Darcy_elevationMax=Darcy_elevation; % i.e. maximum reservoir cond. as soon as soil level underneath reservoir level
% Darcy_elevationMin=[ nan nan PARA.ensemble.initial_altitude(1) ] ; % i.e. altitude of first tile sets

PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1 || strcmp(boundaryCondition{i},'DarcyReservoirNoInflow')==1
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
%     elseif strcmp(boundaryCondition{i},'DarcyReservoirDynamicConductivity')==1
%         PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
%         PARA.ensemble.boundaryCondition(i).parameters.fluxFactorMax=Darcy_fluxFactorMax(i); 
%         PARA.ensemble.boundaryCondition(i).parameters.fluxFactorMin=Darcy_fluxFactorMin(i); 
%         PARA.ensemble.boundaryCondition(i).parameters.elevationMax=Darcy_elevationMax(i);  
%         PARA.ensemble.boundaryCondition(i).parameters.elevationMin=Darcy_elevationMin(i);        
    end
end


% parameters related to snow exchange
% to be specificed by user
%PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);
PARA.ensemble.immobile_snow_height = [ 0.1, 0.1, 0.1 ];
PARA.ensemble.snow_scaling = ones(1, numlabs);  % unclear if needed in ensemble struct


% parameters related to infiltration scheme
% to be specificed by user
% PARA.ensemble.rootDepth = [0.2, 0.2, 0.2 ] ;%[0.2, 0.1, 0.2 ] ;
% PARA.ensemble.fieldCapacity = SETUP.fieldCapacity .* ones(1, numlabs);
% PARA.ensemble.external_water_flux = zeros(1, numlabs);

% location-specific fix parameter values
PARA.location.initial_altitude = PARA.ensemble.initial_altitude(index);
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

% different stratigraphies
% to be specificed by user
PARA.soil.layer_properties = SETUP.stratigraphy;
PARA.soil.layer_properties = PARA.soil.layer_properties{index};

% different initial conditions
% to be specificed by user
PARA.Tinitial = SETUP.Tinitial{labindex};


end