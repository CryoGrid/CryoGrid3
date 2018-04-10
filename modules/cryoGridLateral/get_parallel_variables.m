function PARA = get_parallel_variables(PARA)

index = labindex;

% geometric relations
PARA.ensemble.distanceBetweenPoints=10 .* ( diag(ones(numlabs-1,1),-1)+diag(ones(numlabs-1,1),1) ); %   %in m. Put 0 for all non-connected ensemble members
PARA.ensemble.weight = [1, 1];
PARA.ensemble.area = PARA.ensemble.weight.*100; % in m^2

% topographical relations
PARA.ensemble.initial_altitude = [20.0, 20.5];	%in m a.s.l., this is the reference for the "zero" position of the grids
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;

% parameters related to heat exchange
PARA.ensemble.thermal_contact_length = 20 .* diag(1*ones(numlabs-1,1),-1)+diag(1*ones(numlabs-1,1),1);
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;

% parameters related to water exchange
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per output interval from each worker to worker index
PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * ( ones(numlabs) - eye(numlabs ) );	%in m/sec % [Roth: 1e-5 for saturated silt, 		2.2e-5 for saturated sand]
% Leos definition:  PARA.ensemble.hydraulic_conductivity= diag(PARA.soil.hydraulic_conductivity*ones(numlabs-1,1),-1)+diag(PARA.soil.hydraulic_conductivity*ones(numlabs-1,1),1); %in m/sec % [Roth: 1e-5 for saturated silt, 2.2e-5 for saturated sand] % Léo: 10-5 m/s good for peat also
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.hydraulic_contact_length = 20 .* diag( 1 * ones(numlabs-1,1),-1) + diag( 1 * ones(numlabs-1,1),1);
PARA.ensemble.infiltration_altitude = nan(1, numlabs);
PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;

% JAN: not necessary if absoulte limit realization-specific
%PARA.ensemble.infiltration_limit_altitude=min(PARA.ensemble.altitude-PARA.soil.infiltration_limit_depth); % Léo : We should decide if it is a common value for all or if it can varies from a worker to another
%PARA.ensemble.bottomBucketSoilcTIndex = ones(numlabs,1);

boundaryCondition={'DarcyReservoir','NoBC'};
Darcy_elevation=[19.5 NaN ]; % Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers withour this boundary condition
Darcy_fluxFactor=[5*1e-5/50 NaN ]; % Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance    Defined for now like this, lets see if we wantto define it differently
PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1;
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
    end
end

% parameters related to snow exchange
%PARA.ensemble.snow_fluxes = zeros( 1, numlabs );            % total snow flux in [m SWE] per output interval from each worker to worker index
PARA.ensemble.immobile_snow_height = [0.1, 0.1 ]; %in m
PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);

% parameters related to infiltration scheme
PARA.ensemble.external_water_flux=[0, 0] ;  %in m/day
PARA.ensemble.rootDepth = [ 0.2, 0.1 ];
PARA.ensemble.fieldCapacity = [ 0.5, 0.5 ];

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
% location-specific dynamic common thresholds
PARA.location.absolute_maxWater_altitude = [max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater];
PARA.location.absolute_maxSnow_altitude = [max( PARA.ensemble.altitude ) + PARA.snow.relative_maxSnow];

% JAN: not necessary if absoulte limit realization-specific
%PARA.location.bottomBucketSoilcTIndex = PARA.ensemble.bottomBucketSoilcTIndex(labindex);
%PARA.soil.infiltration_limit_altitude=PARA.ensemble.infiltration_limit_altitude;
    
% different stratigraphies
% to be specificed by user

% different initial conditions
% to be specificed by user
end
