function PARA = get_parallel_variables(PARA)

	index = labindex;

    % auxiliary calculations for circular geometry
    diameter=10; % in [m]
    area = pi.*(diameter./2)^2; % in [m^2]
    perimeter = 2*pi.*(diameter./2); % in [m]

    % geometric relations
    PARA.ensemble.distanceBetweenPoints= diameter .* ( ones(numlabs) - eye(numlabs) );% [0, 4, 4; 4, 0, 4 ; 4, 4, 0];  %in m; put 0 for all non-connected ensemble members
    PARA.ensemble.weight = [1, 1];%[2, 1, 1];  
    PARA.ensemble.area = PARA.ensemble.weight.*area; % in m^2

    % topographical relations
    PARA.ensemble.initial_altitude = [20.0, 20.5]; %[20.0, 21.0, 20.5];                            %in m a.s.l., this is the reference for the "zero" position of the grids
    PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
    PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;

    % parameters related to heat exchange
    PARA.ensemble.thermal_contact_length = perimeter .* ( ones(numlabs) - eye(numlabs ) ); % [ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %
    
    % parameters related to water exchange
    PARA.ensemble.water_fluxes = zeros( 1, numlabs ); % total water flux in [m] per output interval from each worker to worker index
    PARA.ensemble.external_water_flux=[0, 0 ] ; % 0];   %in m/day
    PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * ( ones(numlabs) - eye(numlabs ) );%[ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %in m/sec % [Roth: 1e-5 for saturated silt, 		2.2e-5 for saturated sand]
    PARA.ensemble.water_table_altitude = PARA.ensemble.altitude;  %initialize somehow;    

    %PARA.ensemble.max_water_flux= [0 0];   %in m water equivalent
    PARA.ensemble.hydraulic_contact_length = perimeter .*  ( ones(numlabs) - eye(numlabs ) );%[ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ];
    PARA.ensemble.active_layer_depth_altitude = [NaN NaN NaN];
    PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;
    
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
    
    % different stratigraphies
    PARA.soil.layer_properties = {[0.0    0.5    0.5    0.00   1   0.50 ;...            % center stratigraphy without excess ice
                                   1.0    0.5    0.5    0.00   1   0.50 ;...
                                  10.0    0.25   0.75   0.00   1   0.25     ] , ...
                                  [0.0     0.5    0.5     0.00   1   0.50;...          % rim stratigraphy with excess ice
                                   0.7     0.8    0.2     0.00   1   0.50;...
                                  10.0    0.25   0.75    0.00   1   0.25     ], ...
                                  [0.0     0.5    0.5     0.00   1   0.50;...          % trough stratigraphy with excess ice
                                   0.1     0.8    0.2     0.00   1   0.50;...
                                  10.0    0.25   0.75    0.00   1   0.25     ]};
                              
    PARA.soil.layer_properties = PARA.soil.layer_properties{index};
    
 
    
    % different initial conditions
     PARA.Tinitial = [-5     5    5     5;...
                       0    -5   -5     -5;...
                       1    -5   -5     -5;...
                      10    -8   -8     -8;...
                      20   -10  -10     -10;...
                     100   -10  -10     -10;...  
                    2000    10   10     10];

     PARA.Tinitial=[PARA.Tinitial(:,1) PARA.Tinitial(:, 1+index)];

end
