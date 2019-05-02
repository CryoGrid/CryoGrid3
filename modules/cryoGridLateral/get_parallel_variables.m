% in this function all parameters which are different for the different
% tiles, are specified and overwrite previously set values; in addtion, the
% topological relations are specified
function PARA = get_parallel_variables(PARA)

index = labindex;

% Load geometrical setup
geomSetup = iLoadGeomSetup( 'geomSetup.mat', PARA.ensemble.geomSetup );

% topological relations
PARA.ensemble.weight = geomSetup.weight;
PARA.ensemble.area = geomSetup.area; % in m^2

% For in-line workers
PARA.ensemble.distanceBetweenPoints= geomSetup.distanceBetweenPoints; %   %in m. Put 0 for all non-connected ensemble members
A = geomSetup.A; % adjacency matrix of the network (auxiliary)

% topographical relations
PARA.ensemble.initial_altitude = geomSetup.initial_altitude;	%in m a.s.l., this is the reference for the "zero" position of the grids
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to HEAT exchange
PARA.ensemble.thermal_contact_length = geomSetup.thermal_contact_length;
PARA.ensemble.thermalDistance = geomSetup.thermalDistance;

% parameters related to WATER exchange
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * A;
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
PARA.ensemble.infiltration_altitude = nan(1, numlabs);
PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;

boundaryCondition=geomSetup.boundaryCondition; 		% set to 'DarcyReservoir' for an external water reservoir
Darcy_elevation=geomSetup.Darcy_elevation; 				% Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers without this boundary condition
Darcy_fluxFactor=geomSetup.Darcy_fluxFactor; 			% Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance. Defined for now like this, lets see if we wantto define it differently. NaN for workers without this boundary condition
PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i);
    end
end

% parameters related to snow exchange
% to be specificed by user
PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);
PARA.ensemble.immobile_snow_height = geomSetup.immobile_snow_height;
PARA.ensemble.snow_scaling = ones(1, numlabs);  % New variable for the new snow scheme from Jan

% parameters related to infiltration scheme
% to be specified by user

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


% Stratigraphy and initial Temperature
precalc=1; % Use this guy like a switch depending on if you want to use precalc values or the old school way to calc these values

if precalc==1
    assert(~isempty(geomSetup.thermalInit),'geomSetup.thermalInit is empty')
    PARA.soil.layer_properties=geomSetup.thermalInit(labindex).layer_properties;
    PARA.Tinitial=geomSetup.thermalInit(labindex).Tinitial;
else
    % Active layers
    ActiveLayer=[NaN 0.9 0.9 0.85 0.8 0.75 0.7]; % Input active layers
    ispf=~isnan(ActiveLayer); % Define who is permafrost
    
    % Stratigraphy
    Strati_mire =[    0.0     0.80    0.05    0.15    1   0.80    ;...
                      0.5     0.80    0.05    0.15    1   0.80    ;...
                      3.0     0.50    0.50    0.00    2   0.50    ;...
                     10.0     0.03    0.97    0.00    1   0.03   ];
    Strati_palsa_initial=Strati_mire; % Start from the mire strati
    Strati_palsa_initial(1,2)=PARA.soil.fieldCapacity; % Set the upper layer at FC
    palsaHeight=PARA.ensemble.initial_altitude-min(PARA.ensemble.initial_altitude);
    if ispf(labindex)==0;
        PARA.soil.layer_properties=Strati_mire;
    else
        PARA.soil.layer_properties=stratiXice(Strati_palsa_initial, palsaHeight(labindex), ActiveLayer(labindex)); % Modify the ice, mineral and organic content
    end
    
    % Initial T profiles
    T_z    =[-5 0 0.1 0.5 1 2 10 30 500 1000]';
    ActiveLayer(isnan(ActiveLayer))=0.5;
    if min(ActiveLayer)>0.1 && max(ActiveLayer)<1;
        T_z(4)=ActiveLayer(labindex);
    else
        error('Check initial active layers and initial T profile')
    end
    T_mire =[10 5 2 1.5 2 2 2 2 4 10]';
    T_palsa=[10 5 2 0 -1 -1 1 2 4 10]';
    if ispf(labindex)==0;
        PARA.Tinitial=[T_z T_mire];
    else
        PARA.Tinitial=[T_z T_palsa];
    end
end

end