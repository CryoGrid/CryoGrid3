function PARA = get_parallel_variables(PARA, SETUP)

index = labindex;

% geometric relations ( we assume a 2D geometry ( x - lateral, z - vertical and infinite extent in y-direction )

% set the relative lateral extents of the tiles (as INTEGERS!)
w_plateau = 100;
w_intermediate = 10;
w_cliff = 1;

D = 1;              % use this to scale the lateral extents (x-direction)

% topography
PARA.ensemble.initial_altitude = [ 20.0, 20.0, 20.0 ];  %in m a.s.l., this is the reference for the "zero" position of the grids

% external reservoir
K_Reservoir = 5e-5;
boundaryCondition={'NoBC','NoBC', 'DarcyReservoir'};
elevation_Reservoir = [ nan nan -20.0];                 % specify the elevation of a reservoir relative to the tile's initial altitude ( put nan if no connection to external reservoir, should be consistent with "boundaryCondition")


% from here on everything is calculated based on the values above
L = 1;                                                  % scaling of the y-direction (this is just a dummy which cancels out during further calculations of lateral fluxes)

% areas and weights
PARA.ensemble.weight = [ w_plateau, w_intermediate, w_cliff ];
PARA.ensemble.area = PARA.ensemble.weight .* D .* L;

% thermal distances
PARA.ensemble.distanceBetweenPoints = D.* [ 0, (w_plateau+w_intermediate)./2, 0; (w_plateau+w_intermediate)./2, 0, (w_intermediate+w_cliff)./2; 0, (w_intermediate+w_cliff)./2, 0 ];
A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)

% topographical relations
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to HEAT exchange
PARA.ensemble.thermal_contact_length = L .* A;
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;
PARA.ensemble.thermal_contact_DepthOverDistance = 0; % this factor is used to control the depth below which heat exchange takes place: contact_altitude = min(altitudes) - DepthOverDistance * Distance_ij

% parameters related to WATER exchange
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.infiltration_altitude = nan(1, numlabs);
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
PARA.ensemble.hydraulic_conductivity = PARA.soil.hydraulic_conductivity * A;
PARA.ensemble.hydraulic_contact_length = L .* A;
PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;

Darcy_elevation = PARA.ensemble.initial_altitude + elevation_Reservoir; % Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers withour this boundary condition
Darcy_fluxFactor= ~isnan(Darcy_elevation) .* K_Reservoir; % Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance    Defined for now like this, lets see if we wantto define it differently
PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1;
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
    end
end

% parameters related to SNOW exchange
PARA.ensemble.immobile_snow_height = [ 0.1, 0.1, 0.1 ]; %in m
PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);


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
stratigraphyMap = containers.Map( { 'PLATEAU', 'INTERMEDIATE', 'CLIFF' }, ...
    { [ 0.0   0.30    0.10    0.10    1   0.80;...
        0.15  0.70    0.25    0.05    1   0.70;...
        0.55  0.90    0.05    0.05    1   0.40;...
        20.0  0.30    0.70    0.00    1   0.30     ], ...
    [   0.0   0.30    0.10    0.10    1   0.80;...
        0.15  0.70    0.25    0.05    1   0.70;...
        0.55  0.90    0.05    0.05    1   0.40;...
        20.0  0.30    0.70    0.00    1   0.30     ], ...
    [   0.0   0.30    0.10    0.10    1   0.80;...
        0.15  0.70    0.25    0.05    1   0.70;...
        0.55  0.90    0.05    0.05    1   0.40;...
        20.0  0.30    0.70    0.00    1   0.30     ], ...
    } );

PARA.soil.layer_properties = { stratigraphyMap('PLATEAU'), ...
    stratigraphyMap('INTERMEDIATE'), ...
    stratigraphyMap('CLIFF') };

PARA.soil.layer_properties = PARA.soil.layer_properties{index};

end
