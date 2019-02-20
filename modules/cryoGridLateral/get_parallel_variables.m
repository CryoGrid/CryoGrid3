% in this function all parameters which are different for the different
% tiles, are specified and overwrite previously set values; in addtion, the
% topological relations are specified
function PARA = get_parallel_variables(PARA)

index = labindex;

% topological relations
area_tot = 200 + 40 + 40 + 40 + 200;
PARA.ensemble.weight = [20 4 4 4 20]; % Has to be an integer
PARA.ensemble.area = PARA.ensemble.weight ./ sum(PARA.ensemble.weight) .* area_tot ; % in m^2
% For in-line workers
dist=[1.20 0.20 0.20 1.20];
A=zeros(numlabs,numlabs);
idx = sub2ind(size(A),[1:numlabs-1 2:numlabs],[2:numlabs 1:numlabs-1]);
A(idx) = [dist dist];
PARA.ensemble.distanceBetweenPoints= A; %   %in m. Put 0 for all non-connected ensemble members
A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)

% topographical relations
PARA.ensemble.initial_altitude = [300.0000  300.2500  300.5000  300.7500  301.0000];	%in m a.s.l., this is the reference for the "zero" position of the grids
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to HEAT exchange
PARA.ensemble.thermal_contact_length = 5 .* A;
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;

% parameters related to WATER exchange
PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
PARA.ensemble.hydraulic_conductivity= PARA.soil.hydraulic_conductivity * A;
PARA.ensemble.water_table_altitude = nan(1, numlabs);
PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
PARA.ensemble.infiltration_altitude = nan(1, numlabs);
PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;

boundaryCondition={'DarcyReservoir','NoBC','NoBC','NoBC','NoBC'}; 		% set to 'DarcyReservoir' for an external water reservoir
Darcy_elevation=[300 nan nan nan nan]; 				% Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers without this boundary condition
Darcy_fluxFactor=[1000*5*1e-5/1.20 nan nan nan nan]; 			% Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance. Defined for now like this, lets see if we wantto define it differently. NaN for workers without this boundary condition
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
PARA.ensemble.immobile_snow_height = [ 0.05 0.05 0.05 0.05 0.05 ];
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


% different stratigraphies
% to be specificed by user
% Start from mire strati
Strati_mire =[    0.0     0.80    0.05    0.15    1   0.80    ;...
                  0.5     0.80    0.05    0.15    1   0.80    ;...
                  3.0     0.50    0.50    0.00    2   0.50    ;...
                 10.0     0.03    0.97    0.00    1   0.03   ];

% Construct Palsa strati for each worker
Strati_palsa_initial=Strati_mire; % Start from the mire strati
Strati_palsa_initial(1,2)=PARA.soil.fieldCapacity; % Set the upper layer at FC
palsaHeight=PARA.ensemble.initial_altitude(2:end)-PARA.ensemble.initial_altitude(1); % Base palsa height on elevations
ActiveLayer=[0.9 0.9 0.8 0.7]; % Input active layers
if labindex < 2;
    PARA.soil.layer_properties=Strati_mire;
else
    PARA.soil.layer_properties=stratiXice(Strati_palsa_initial, palsaHeight(labindex-1), ActiveLayer(labindex-1)); % Modify the ice, mineral and organic content
end


% different initial conditions
% to be specificed by user
T_z    =[-5 0 0.1 0.5 1 2 10 30 500 1000]';
ActiveLayer=[0.5 ActiveLayer];
if min(ActiveLayer)>0.1 && max(ActiveLayer)<1;
    T_z(4)=ActiveLayer(labindex);
else
    error('Check initial active layers and initial T profile')
end
T_mire =[10 5 2 1.5 2 2 2 2 4 10]';
T_palsa=[10 5 2 0 -1 -1 1 2 4 10]';
if labindex < 2;
    PARA.Tinitial=[T_z T_mire];
else
    PARA.Tinitial=[T_z T_palsa];
end

end