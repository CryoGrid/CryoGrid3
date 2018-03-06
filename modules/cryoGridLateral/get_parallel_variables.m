function PARA = get_parallel_variables(PARA)

	index = labindex;

    % auxiliary calculations for circular geometry
    diameter=10; % in [m]
    area = pi.*(diameter./2)^2; % in [m^2]
    perimeter = 2*pi.*(diameter./2); % in [m]

    % geometric relations
    PARA.ensemble.distanceBetweenPoints=[0 2; 2 0]; %[0, 4, 4; 4, 0, 4 ; 4, 4, 0];% diameter .* ( ones(numlabs) - eye(numlabs) );%   %in m; put 0 for all non-connected ensemble members
    PARA.ensemble.weight = [1, 1];%[2, 1, 1];  
    PARA.ensemble.area = [1 1]; % PARA.ensemble.weight.*area; % in m^2

    % topographical relations
    PARA.ensemble.initial_altitude = [300 301.25]; %[20.0, 21.0, 20.5];                            %in m a.s.l., this is the reference for the "zero" position of the grids
    PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
    PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;

    % parameters related to heat exchange
    PARA.ensemble.thermal_contact_length = perimeter .* ( ones(numlabs) - eye(numlabs ) ); % [ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %
    
    % parameters related to water exchange
    PARA.ensemble.water_fluxes = zeros( 1, numlabs ); % total water flux in [m] per output interval from each worker to worker index
    PARA.ensemble.external_water_flux=[0 0 ] ; % 0];   %in m/day
    PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * ( ones(numlabs) - eye(numlabs ) );%[ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %in m/sec % [Roth: 1e-5 for saturated silt, 		2.2e-5 for saturated sand] % Léo: 10-5 m/s good for peat also
    PARA.ensemble.water_table_altitude = PARA.ensemble.altitude;  %initialize somehow;
    PARA.ensemble.alt_infiltration_limit=min(PARA.ensemble.altitude-PARA.soil.infiltration_limit); % *ones(length(PARA.ensemble.area),1); % LEO: We should decide if it is a common value for all or if it can varies from a worker to another
    %PARA.ensemble.max_water_flux= [0 0];   %in m water equivalent
    PARA.ensemble.hydraulic_contact_length = 1 .*  ( ones(numlabs) - eye(numlabs ) );%[ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ];
    PARA.ensemble.active_layer_depth_altitude = [NaN NaN];
    PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;
    PARA.ensemble.bottomBucketSoilcTIndex = [1 1];
    boundaryCondition={'DarcyReservoir','NoBC'};
    Darcy_elevation=[300 NaN]; % Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers withour this boundary condition
    Darcy_fluxFactor=[5*1e-5/50 NaN]; % Taken as the section hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance    Defined for now like this, lets see if we wantto define it differently
    PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
    [PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
    for i=1:length(boundaryCondition)
        if strcmp(boundaryCondition{i},'DarcyReservoir')==1;
            PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
            PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
        end
    end
    
    
    % parameters related to snow exchange
    %PARA.ensemble.snow_diffusivity = PARA.snow.diffusivity;
    %PARA.ensemble.relative_max_snow_height = 0.2;
    PARA.ensemble.snow_fluxes = zeros( 1, numlabs );            % total snow flux in [m SWE] per output interval from each worker to worker index
    PARA.ensemble.immobile_snow_height = [0.1, 0.1 ]; %, 0.2 ];  %in m %this replaces PARA.snow.maxSnow ?
    PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);
    PARA.ensemble.snow_contact_length = perimeter .* ( ones(numlabs) - eye(numlabs ) );

    % location-specific static
    PARA.location.initial_altitude = PARA.ensemble.initial_altitude(index);
    PARA.soil.externalWaterFlux = PARA.ensemble.external_water_flux(index);
	% location-specific dynamic auxiliary variables
    PARA.location.area = PARA.ensemble.area(index);
    PARA.location.altitude = PARA.ensemble.altitude(index);
    PARA.location.surface_altitude = PARA.ensemble.surface_altitude(index);
    PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(index);
	PARA.location.active_layer_depth_altitude = PARA.ensemble.active_layer_depth_altitude(index);
	% location-specific dynamic common thresholds
	PARA.location.absolute_maxWater_altitude = [max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater];
    PARA.location.absolute_maxSnow_altitude = [max( PARA.ensemble.altitude ) + PARA.snow.relative_maxSnow];
    PARA.location.bottomBucketSoilcTIndex = PARA.ensemble.bottomBucketSoilcTIndex(labindex);
    
    % different stratigraphies
    layer_properties=[ 0.0   0.55    0.05    0.15   1   0.80;...
                       0.5   0.80    0.05    0.15   2   0.80;...
                       3.0   0.50    0.50    0.0    1   0.50;...
                      10.0   0.03    0.97    0.0    1   0.03 ];
    
    PARA.soil.layer_properties = { layer_properties, layer_properties};
                              
    PARA.soil.layer_properties = PARA.soil.layer_properties{index};
    
 
    
    % different initial conditions
PARA.Tinitial = [-5     10   10  ;...
                  0     5    5   ;...
                  0.1   2    2   ;...
                  0.5   1.2  0.5 ;... 
                  1     1    0   ;...  
                  2     0.8 -0.2 ;...
                  10    1    0   ;...
                  30    2    2   ;...
                  500   4    4   ;...
                  5000  10   10];

     PARA.Tinitial=[PARA.Tinitial(:,1) PARA.Tinitial(:, 1+index)];

end
