function PARA = get_parallel_variables(PARA,SETUP)

	index = labindex;
    
    % auxiliary calculations for circular geometry
    L_F = SETUP.LF; % landscape Lake Fraction
    radius = SETUP.LR; % Lake radius in [m]
    perimeter = 2*pi.*radius; % in [m]
    distance =  sqrt(pi./(4*L_F)) * radius; % in [m]
    %distance = 2.*radius;
    % geometric relations
%    PARA.ensemble.distanceBetweenPoints = diameter .* ( ones(numlabs) - eye(numlabs) );% [0, 4, 4; 4, 0, 4 ; 4, 4, 0];  %in m; put 0 for all non-connected ensemble members
    PARA.ensemble.distanceBetweenPoints = distance .* ( ones(numlabs) - eye(numlabs) );  %in m; put 0 for all non-connected ensemble members;     %zzz needed?
    %%%PARA.ensemble.distanceBetweenPoints = distance ;  %in m; put 0 for all non-connected ensemble members;     %zzz needed?
    %    PARA.ensemble.weight = [1, 1];%[2, 1, 1];  
    %%%PARA.ensemble.weight = [radius/(PARA.ensemble.distanceBetweenPoints+radius),PARA.ensemble.distanceBetweenPoints/(PARA.ensemble.distanceBetweenPoints+radius)]; % weight keff by length radius and distance  
    PARA.ensemble.weight = [radius/(distance+radius),distance/(distance+radius)]; % weight keff by length radius and distance  

    %area = pi.*radius^2; % in [m^2]
    %area = [pi.*radius^2 , pi.*(PARA.ensemble.distanceBetweenPoints^2-radius^2); % in [m^2]
    %PARA.ensemble.area = PARA.ensemble.weight.*area; % in m^2  zzz jjj ???

    PARA.ensemble.area = [pi.*radius^2 , pi.*(distance^2-radius^2)]; % in [m^2]

    % topographical relations
    %tsvd PARA.ensemble.initial_altitude = [20.0, 20.5]; %[20.0, 21.0, 20.5];                            %in m a.s.l., this is the reference for the "zero" position of the grids
%    PARA.ensemble.initial_altitude = [20.0, 20.0]; %[20.0, 21.0, 20.5];                            %in m a.s.l., this is the reference for the "zero" position of the grids
    PARA.ensemble.initial_altitude = [0., 0.]; %[20.0, 21.0, 20.5];                            %in m a.s.l., this is the reference for the "zero" position of the grids

    PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
    PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;

    % parameters related to heat exchange
    PARA.ensemble.thermal_contact_length = perimeter .* ( ones(numlabs) - eye(numlabs ) ); % [ 0, 1, 0 ; 1, 0, 1 ; 0, 1, 0 ]; %
    
    % parameters related to water exchange
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
    PARA.ensemble.immobile_snow_height = [0.1, 0.1 ]; %, 0.2 ];  %in m %this replaces PARA.snow.maxSnow ? zzz
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
    
    % soil stratigraphy
    % column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer extends until the end of the model domain
    % column 2: volumetric water+ice content; column 3: volumetric mineral content; column 4: volumetric organic content;
    % column 5: code for soil type: 1: sand, 2: silt
    % column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs   
    
    
    % default stratigraphy used in publication:
PARA.soil.layer_properties={[    0.0   0.60    0.10    0.15    1   0.75;...         % Non-Lake
                                 0.15  0.65    0.3     0.05    2   0.65;...
                                 0.9   0.65    0.3     0.05    1   0.65;...
                                 9.0   0.30    0.70    0.00    1   0.30     ], ...
                                %
                            [    0.0   0.65    0.3     0.05    2   0.65;...          % Lake
                                 0.9   0.65    0.3     0.05    1   0.65;...
                                 9.0   0.30    0.70    0.00    1   0.30     ]};
                                                                
 PARA.soil.layer_properties = PARA.soil.layer_properties{index};
      
    % different initial conditions             

     PARA.Tinitial=[SETUP.Tini(:,1) SETUP.Tini(:,1+index)];
     
     PARA.water.depth = [0.,SETUP.LD];
     PARA.water.depth = PARA.water.depth(index); 
end
