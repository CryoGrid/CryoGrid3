function PARA = get_parallel_variables(PARA)

	index = labindex;
    
    % auxiliary calculations for circular geometry
    F_L = 0.5; % landscape Lake Fraction
    radius = 10; % in [m]
    perimeter = 2*pi.*radius; % in [m]
    distance =  sqrt(pi./(4*F_L)) * radius; % in [m]
    % geometric relations
%    PARA.ensemble.distanceBetweenPoints = diameter .* ( ones(numlabs) - eye(numlabs) );% [0, 4, 4; 4, 0, 4 ; 4, 4, 0];  %in m; put 0 for all non-connected ensemble members
    PARA.ensemble.distanceBetweenPoints = distance .* ( ones(numlabs) - eye(numlabs) );  %in m; put 0 for all non-connected ensemble members;     %zzz needed?
    PARA.ensemble.weight = [1, 1];%[2, 1, 1];  
%    PARA.ensemble.weight = [radius/(PARA.ensemble.distanceBetweenPoints+radius),PARA.ensemble.distanceBetweenPoints/(PARA.ensemble.distanceBetweenPoints+radius)]; % weight keff by length radius and distance  

 %  area = pi.*(diameter./2)^2; % in [m^2]
    %area = [pi.*radius^2 , pi.*(PARA.ensemble.distanceBetweenPoints^2-radius^2); % in [m^2]
%    PARA.ensemble.area = PARA.ensemble.weight.*area; % in m^2  zzz jjj ???
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
%     PARA.soil.layer_properties = {[0.0    0.5    0.5    0.00   1   0.50 ;...            % center stratigraphy without excess ice
%                                    1.0    0.5    0.5    0.00   1   0.50 ;...
%                                   10.0    0.25   0.75   0.00   1   0.25     ] , ...
%                                   [0.0     0.5    0.5     0.00   1   0.50;...          % rim stratigraphy with excess ice
%                                    0.7     0.8    0.2     0.00   1   0.50;...
%                                   10.0    0.25   0.75    0.00   1   0.25     ], ...
%                                   [0.0     0.5    0.5     0.00   1   0.50;...          % trough stratigraphy with excess ice
%                                    0.1     0.8    0.2     0.00   1   0.50;...
%                                   10.0    0.25   0.75    0.00   1   0.25     ]};


  PARA.soil.layer_properties = {[0.0    0.5    0.5    0.00   1   0.50 ;...          % non-lake stratigraphy (no excess ice)
                                 1.0    0.5    0.5    0.00   1   0.50 ;...
                                10.0    0.25   0.75   0.00   1   0.25     ], ...
                                  
                                [0.0    0.5    0.5    0.00   1   0.50 ;...          % lake stratigraphy 
                                 1.0    0.5    0.5    0.00   1   0.50 ;...
                                10.0    0.25   0.75   0.00   1   0.25     ]};
 
%todotodo  read parameters from file!
%  PARA.soil.layer_properties = {[0.0    0.5    0.5    0.25   1   0.50 ;...            % reference case - Langer 2015
%                                 0.2    0.5    0.5    0.00   1   0.50 ;...
%                                 
%                                 ... update
%                                 0.9    0.25   0.75   0.00   1   0.25     ] , ...
%                                 9.0    0.25   0.75   0.00   1   0.25     ] , ...
%                              1000.0    0.25   0.75   0.00   1   0.25     ] , ...
                                  
                              
                              %     PARA.soil.layer_properties = {[0.0     0.5    0.5    0.00   1   0.50 ;...            % lake stratigraphy (no excess ice)
%                                    0.02    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.04    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.06    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.08    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.10    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.12    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.14    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.16    0.5    0.5    0.00   1   0.50 ;...
%                                    0.18    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.20    0.5    0.5    0.00   1   0.50 ] , ...
%                                    
%                                   [0.0     0.5    0.5    0.00   1   0.50 ;...            % lab 2
%                                    0.02    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.04    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.06    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.08    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.10    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.12    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.14    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.16    0.5    0.5    0.00   1   0.50 ;...
%                                    0.18    0.5    0.5    0.00   1   0.50 ;... 
%                                    0.20    0.5    0.5    0.00   1   0.50 ] };                                 
 PARA.soil.layer_properties = PARA.soil.layer_properties{index};
      
    % different initial conditions
     PARA.Tinitial = [-5     5    5;...
                       0    -5   -2;...
                       1    -5   -1;...
                      10    -8    0;...
                      20   -10    0;...
                     100   -10    1;...  
                    2000    10   10];
%     PARA.Tinitial = [-5   10   10;...    this profile can cause model crashes in mpi...
%                       0    0    0;...
%                       5   -5   -5;...
%                       20   -10 -10;...
%                       100  -10 -10;...
%                       2000  10  10];                

     PARA.Tinitial=[PARA.Tinitial(:,1) PARA.Tinitial(:, 1+index)];
     
%tsvd  set lake depths for workers     
%     PARA.water.depth = [1.,0.]; 

%     PARA.water.depth = [0.,0.]; % non-lake setting 
     PARA.water.depth = [0.,0.9]; % non-lake and small-sized water body (SSW) 
     %PARA.water.depth = [0.,6.]; % non-lake and medium-sized water body (MSW) 

     PARA.water.depth = PARA.water.depth(index); 
end
