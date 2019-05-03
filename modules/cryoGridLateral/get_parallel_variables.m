% in this function all parameters which are different for the different
% tiles, are specified and overwrite previously set values; in addtion, the
% topological relations are specified
function PARA = get_parallel_variables( PARA )

index = labindex;

% the following choice of topological parameters is just an example for three connected linear tiles

PARA.ensemble.area = [ 10 10 10 ];
PARA.ensemble.weight = PARA.ensemble.area;

PARA.ensemble.distanceBetweenPoints = [ 0, 1, 0; 1, 0, 1; 0, 1, 0 ];

% parameters related to HEAT exchange
PARA.ensemble.thermal_contact_length = [0, 10, 0; 10, 0, 10; 0, 10, 0 ];
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;

% parameters related to WATER exchange
PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
PARA.ensemble.hydraulicDistance = PARA.ensemble.thermal_contact_length;

% topographical relations
A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)

PARA.ensemble.initial_altitude = [ 0, 1, 2 ]; %in m a.s.l.
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

boundaryCondition={'NoBC', 'NoBC', 'NoBC'}; 		% set this to 'DarcyReservoir' for an external water reservoir
Darcy_elevation=[ nan, nan, nan ] ; 				% Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers without this boundary condition
Darcy_fluxFactor= [ nan nan nan ];					% Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance. Defined for now like this, lets see if we wantto define it differently. NaN for workers without this boundary condition


PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1 || strcmp(boundaryCondition{i},'DarcyReservoirNoInflow')==1
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i);        
    end
end

% parameters related to snow exchange
% to be specificed by user
PARA.ensemble.immobile_snow_height = [ 0.1, 0.1, 0.1 ];
PARA.ensemble.snow_scaling = ones(1, numlabs);

% parameters related to infiltration scheme
% to be specificed by user
% ...

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
% ...

% different initial conditions
% to be specificed by user
% ...
end
