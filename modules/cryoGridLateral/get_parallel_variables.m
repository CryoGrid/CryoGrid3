function PARA = get_parallel_variables(PARA)
% Case : The long slope

    % auxiliary calculations for circular geometry
    diameter=10; % in [m]
    area = pi.*(diameter./2)^2; % in [m^2]
    perimeter = 4; % 2*pi.*(diameter./2); % in [m]

    % geometric relations
    PARA.ensemble.distanceBetweenPoints=diag(2*ones(numlabs-1,1),-1)+diag(25*ones(numlabs-1,1),1); %   %in m. Put 0 for all non-connected ensemble members
    PARA.ensemble.weight = [1 1 1 1 1];
    PARA.ensemble.area = [1 1 1 1 1]; % in m^2

    % topographical relations
    PARA.ensemble.initial_altitude = [300 300.3 300.6 300.9 301.2]; %in m a.s.l., this is the reference for the "zero" position of the grids
    PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
    PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;

    % parameters related to heat exchange
    PARA.ensemble.thermal_contact_length = diag(1*ones(numlabs-1,1),-1)+diag(1*ones(numlabs-1,1),1);
    
    % parameters related to water exchange
    PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per output interval from each worker to worker index
    PARA.ensemble.external_water_flux=[0 0 0 0 0] ; % 0];   %in m/day
    PARA.ensemble.hydraulic_conductivity= diag(PARA.soil.hydraulic_conductivity*ones(numlabs-1,1),-1)+diag(PARA.soil.hydraulic_conductivity*ones(numlabs-1,1),1); %in m/sec % [Roth: 1e-5 for saturated silt, 2.2e-5 for saturated sand] % Léo: 10-5 m/s good for peat also
    PARA.ensemble.water_table_altitude = PARA.ensemble.altitude;  %initialize somehow;
    PARA.ensemble.infiltration_limit_altitude=min(PARA.ensemble.altitude-PARA.soil.infiltration_limit_depth); % Léo : We should decide if it is a common value for all or if it can varies from a worker to another
    %PARA.ensemble.max_water_flux= [0 0];   %in m water equivalent
    PARA.ensemble.hydraulic_contact_length = diag( 1 * ones(numlabs-1,1),-1) + diag( 1 * ones(numlabs-1,1),1);
    PARA.ensemble.infiltration_altitude = nan(numlabs,1);
    PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;
    PARA.ensemble.bottomBucketSoilcTIndex = ones(numlabs,1);
    boundaryCondition={'DarcyReservoir','NoBC', 'NoBC','NoBC','NoBC'};
    Darcy_elevation=[300 NaN NaN NaN NaN]; % Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers withour this boundary condition
    Darcy_fluxFactor=[5*1e-5/50 NaN NaN NaN NaN]; % Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance    Defined for now like this, lets see if we wantto define it differently
    PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
    [PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
    for i=1:numlabs
        if strcmp(boundaryCondition{i},'DarcyReservoir')==1;
            PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
            PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
        end
    end
    
    
    % parameters related to snow exchange
    %PARA.ensemble.snow_diffusivity = PARA.snow.diffusivity;
    %PARA.ensemble.relative_max_snow_height = 0.2;
    PARA.ensemble.snow_fluxes = zeros( 1, numlabs );            % total snow flux in [m SWE] per output interval from each worker to worker index
    PARA.ensemble.immobile_snow_height = [0.1, 0.1, 0.1, 0.1, 0,1]; %, 0.2 ];  %in m %this replaces PARA.snow.maxSnow ?
    PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);
    PARA.ensemble.snow_contact_length = PARA.ensemble.thermal_contact_length; % perimeter .* ( ones(numlabs) - eye(numlabs ) );

    % location-specific static
    PARA.location.initial_altitude = PARA.ensemble.initial_altitude(labindex);
    PARA.soil.externalWaterFlux = PARA.ensemble.external_water_flux(labindex);
	% location-specific dynamic auxiliary variables
    PARA.location.area = PARA.ensemble.area(labindex);
    PARA.location.altitude = PARA.ensemble.altitude(labindex);
    PARA.location.surface_altitude = PARA.ensemble.surface_altitude(labindex);
    PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
	PARA.location.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
	% location-specific dynamic common thresholds
	PARA.location.absolute_maxWater_altitude = max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater; % Leo : change max->min if you prefer that a relative value of 0, block ponding for everybody. 
    PARA.location.absolute_maxSnow_altitude = max(PARA.ensemble.altitude) + PARA.snow.relative_maxSnow;
    PARA.location.bottomBucketSoilcTIndex = PARA.ensemble.bottomBucketSoilcTIndex(labindex);
    PARA.soil.infiltration_limit_altitude=PARA.ensemble.infiltration_limit_altitude;
    
    % different stratigraphies
    layer_properties=[ 0.0   0.55    0.05    0.15   1   0.80;...
                       0.5   0.80    0.05    0.15   2   0.80;...
                       3.0   0.50    0.50    0.0    1   0.50;...
                      10.0   0.03    0.97    0.0    1   0.03 ];
    
    PARA.soil.layer_properties = { layer_properties, layer_properties, layer_properties, layer_properties, layer_properties };
                              
    PARA.soil.layer_properties = PARA.soil.layer_properties{labindex};
    
 
    
    % different initial conditions
PARA.Tinitial = [-5     10    10    10    10    10;...
                  0     5     5     5     5     5;...
                  0.1   2     0.5   0.5   0.5   0.5;...
                  0.5   1.2  -1    -1    -1    -1;... 
                  1     1    -1    -1    -1    -1;...  
                  2     0.8  -0.5  -0.5  -0.5  -0.5;...
                  10    1     0     0     0     0;...
                  30    2     2     2     2     2;...
                  500   4     4     4     4     4;...
                  5000  10    10    10    10    10 ];
              
PARA.Tinitial=[PARA.Tinitial(:,1) PARA.Tinitial(:, 1+labindex)];

end
