% in this function all parameters which are different for the different
% tiles, are specified and overwrite previously set values; in addtion, the
% topological relations are specified
function PARA = get_parallel_variables( PARA, SETUP )

[mesoIndex, microIndex] = getMesoMicroIndices(labindex, PARA.ensemble.microUnitSizes);

PARA.soil.hydraulic_conductivity_subs = SETUP.Ksubs;
PARA.soil.hydraulic_conductivity_surf_micro = SETUP.Ksurf_micro;
PARA.soil.hydraulic_conductivity_surf_meso = SETUP.Ksurf_meso;
PARA.technical.mesoSurfaceWaterScheme = SETUP.mesoSurfaceWaterScheme;
% parameters related to WATER exchange
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
PARA.ensemble.hydraulic_conductivity_subs = PARA.soil.hydraulic_conductivity_subs .* (PARA.ensemble.adjacency_micro + PARA.ensemble.adjacency_meso);
PARA.ensemble.hydraulic_conductivity_surf = PARA.soil.hydraulic_conductivity_surf_micro .* PARA.ensemble.adjacency_micro + PARA.soil.hydraulic_conductivity_surf_meso .* PARA.ensemble.adjacency_meso;
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.infiltration_altitude = nan(1, numlabs);

boundaryCondition=cell(1,numlabs); 		% set this to 'DarcyReservoir' for an external water reservoir
boundaryCondition(:)={'NoBC'};
Darcy_elevation=nan(1,numlabs) ; 		% Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers without this boundary condition
Darcy_fluxFactor=nan(1,numlabs);		% Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance. Defined for now like this, lets see if we wantto define it differently. NaN for workers without this boundary condition

% boundary condition of first meso unit
boundaryCondition{ PARA.ensemble.microUnitSizes(1) }=SETUP.boundaryConditionFirstMeso;
Darcy_elevation( PARA.ensemble.microUnitSizes(1) )=PARA.ensemble.initial_altitude(1)+SETUP.eReservoirFirstMeso;
Darcy_fluxFactor( PARA.ensemble.microUnitSizes(1) )=2*pi*PARA.soil.hydraulic_conductivity;

% boundary condition of last meso unit
if SETUP.Nmeso>1
    boundaryCondition{ numlabs }= SETUP.boundaryConditionLastMeso;
    Darcy_elevation( numlabs )=PARA.ensemble.initial_altitude( sum( PARA.ensemble.microUnitSizes(1:end-1) )+1  )+SETUP.eReservoirLastMeso;
    Darcy_fluxFactor( numlabs )=2*pi*PARA.soil.hydraulic_conductivity;
end

PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end}; % this is necessary to initialize the right length
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1 || strcmp(boundaryCondition{i},'DarcyReservoirNoInflow')==1
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i);        
    end
end

% parameters related to snow exchange
% to be specificed by user
PARA.ensemble.immobile_snow_height = PARA.snow.catchHeight .* ones(1, numlabs);
PARA.ensemble.snow_scaling = ones(1, numlabs);

% parameters related to infiltration scheme
% to be specificed by user
% ...

% location-specific fix parameter values
PARA.location.initial_altitude = PARA.ensemble.initial_altitude(labindex);
% location-specific dynamic auxiliary variables
PARA.location.area = PARA.ensemble.area(labindex);
PARA.location.altitude = PARA.ensemble.altitude(labindex);
PARA.location.surface_altitude = PARA.ensemble.surface_altitude(labindex);
PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
PARA.location.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
PARA.location.soil_altitude = PARA.ensemble.soil_altitude(labindex);
% location-specific dynamic common thresholds
PARA.location.absolute_maxWater_altitude = max( PARA.ensemble.altitude ) + PARA.soil.relative_maxWater;
PARA.location.absolute_maxSnow_altitude = max( PARA.ensemble.altitude ) + PARA.snow.relative_maxSnow;

% different stratigraphies
% to be specificed by user
PARA.soil.layer_properties = SETUP.stratigraphy{microIndex};

% different initial conditions
% to be specificed by user
PARA.Tinitial = SETUP.Tinitial{labindex};