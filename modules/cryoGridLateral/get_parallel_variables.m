function PARA = get_parallel_variables(PARA)
% in this function all parameters which are different for the different tiles, are specified and overwrite previously set values; 
% in addition, the topological relations are specified

%tsvd  area definition modified!!!  zzz check need to update lateral water exchange etc.
%tsvd  lateral fluxes are only allowed between adjacent workers!

index = labindex; % index 1: tundra, index 2: gravel road

%rel_maxWater=PARA.soil.relative_maxWater; rel_maxWater_pond=PARA.soil.relative_maxWater_ponding;
%tsvd IS area_tot = 70.0;

%tsvd IS  Embankent specifications
EBHag = PARA.IS.EBHag;   % embankment height above ground
EBHbg = PARA.IS.EBHbg;   % embankment height below ground
EBT  =  EBHag+EBHbg;  % total embankment thickness
 
switch numlabs % define tile widths
    case 1
        disp('single mode - 1 tile');
        PARA.ensemble.TileType={'tundra'}; 
        %PARA.ensemble.TileType={'road'}; 
        %PARA.ensemble.TileType={'pile'}; 
        %PARA.ensemble.TileType={'foundation_base'};
        PARA.ensemble.TileWidth = 100;
        disp('single mode - tundra reference');
        PARA.IS.TileType={'tundra'};  
        %PARA.ensemble.TileWidth = 10;
    case 2
        disp('1 gravel road tile, 1 tundra tile');
        PARA.ensemble.TileType = {'road','tundra'};  
        PARA.ensemble.TileWidth = [10 50]; % cccc 7.5 does not work....!
    case 5
        %disp('1 gravel road surface tile, 1 above-ground shoulder tile, 1 below-ground shoulder tile, 1 tundra tile, 1 buffer tundra tile (large extent)');
        %disp('1 gravel road surface tile, 2 above-ground shoulder tile, 1 tundra tile, 1 buffer tundra tile (large extent)');
        disp('2 gravel road surface tiles, 1 shoulder tile, 1 toe tile, 1 tundra tile');
        PARA.ensemble.TileType = {'road','road','shoulder','toe','tundra'};         
        %PARA.ensemble.TileWidth = [5 5 5 10 75];  % tile1: gravel road surface, tile 2: above-ground shoulder, tile 3: below-ground shoulder (toe), tile 4: tundra, tile 5: tundra (large extent)
        PARA.ensemble.TileWidth = [4 1 5 10 40]; % T1: road center, T2: outer road, T3: shoulder, T4: Toe, T5: Tundra (covering 50m distance from the road)
    case {30,36}
        disp([num2str(numlabs),' tiles']); % tile1: gravel road surface, tile 2: above-ground shoulder, tile 3: below-ground shoulder (toe), tile 4: tundra, tile 5: tundra (large extent) 
        PARA.ensemble.TileType(1:5)={'road'};  PARA.ensemble.TileType(6:10)={'shoulder'};  PARA.ensemble.TileType(11:20)={'toe'};  
        PARA.ensemble.TileType(21:numlabs)={'tundra'}; % slope 1:2 (5m shoulder, Embankment height ag 2.5m)
        PARA.ensemble.TileWidth(1:25) =  1.; % 1m tile width for distance 0-25m
        PARA.ensemble.TileWidth(26:27) = 5.;        % 5m tile width for distance 25-35m
        PARA.ensemble.TileWidth(28:29) = 10.;       % 10m distance 35-55m
        PARA.ensemble.TileWidth(30) =    45.;          % 45m distance 55-100m
       % if(numlabs==36); PARA.ensemble.TileWidth(31:36)=...; end 
   case 7 % Diesel Tank Norilsk
        disp('6 Tank tiles, 1 disturbed tundra tile');
        %PARA.ensemble.TileType = {'pile','tank_bottom','tank_bottom','shoulder','toe','foundation_base','disturbed tundra'};         
        PARA.ensemble.TileType = {'tank_bottom','tank_bottom','shoulder','shoulder','toe','foundation_base','disturbed tundra'};         
        PARA.ensemble.TileWidth = [22.5 2.5 2.5 2.5 10 20 40]; % W1: distance from Tank center (Tank diameter: 50m), W2: outer rim (1m), W3: upper shoulder, W4: lower shoulder, W5: tetention basin, W6: foundation base, W7: disturbed tundra 
end

PARA.IS.NumberTiles_Road = sum((string( PARA.ensemble.TileType)=='road')); % number tiles Road Surface
PARA.IS.NumberTiles_Shoulder = sum((string( PARA.ensemble.TileType)=='shoulder')); % number tiles Embankment Shoulder
PARA.IS.NumberTiles_IS = PARA.IS.NumberTiles_Road+PARA.IS.NumberTiles_Shoulder;
PARA.IS.NumberTiles_Toe = sum((string( PARA.ensemble.TileType)=='toe')); % number tiles Toe
PARA.IS.NumberTiles_Tundra = sum((string( PARA.ensemble.TileType)=='disturbed tundra')); % number tiles Tundra

%tsvd IS  PARA.ensemble.distanceBetweenPoints= 10 .* ( diag(ones(numlabs-1,1),-1)+diag(ones(numlabs-1,1),1) ); %   %in m. Put 0 for all non-connected ensemble members
rel_maxSnow = PARA.snow.relative_maxSnow;        
TileDistance = 0.5*(PARA.ensemble.TileWidth(1:length(PARA.ensemble.TileWidth)-1)+PARA.ensemble.TileWidth(2:length(PARA.ensemble.TileWidth))); % distance between mid-points of tiles
DistanceRC = cumsum(PARA.ensemble.TileWidth)-0.5*PARA.ensemble.TileWidth; % distance of tile mid-points to Road Centre
WidthRoad = sum(PARA.ensemble.TileWidth(1:PARA.IS.NumberTiles_Road)); % Road
WidthShoulder = sum(PARA.ensemble.TileWidth(PARA.IS.NumberTiles_Road+1:PARA.IS.NumberTiles_Road+PARA.IS.NumberTiles_Shoulder));
WidthToe = sum(PARA.ensemble.TileWidth(PARA.IS.NumberTiles_Road+PARA.IS.NumberTiles_Shoulder+1:PARA.IS.NumberTiles_Road+PARA.IS.NumberTiles_Shoulder+PARA.IS.NumberTiles_Toe));
SlopeEB = EBHag/WidthShoulder; % embankment slope
SlopeSnow = (EBHag-rel_maxSnow)/(WidthShoulder+WidthToe);

for n=1:numlabs
%% Gravel Road
    if strcmp(PARA.Exp.Case,'GravelRoad')
        switch string(PARA.ensemble.TileType(n))            
            case 'road'
                PARA.ensemble.initial_altitude(n) = EBHag;
                PARA.ensemble.snow.relative_maxSnow(n) = 0.; % max value relative_maxSnow: 0.93m for tile 3 
                PARA.ensemble.soil.albedo(n)=0.3;   % albedo snow-free surface    Andersland, Landanyi  table 3.2  absorptivity of concrete 0.6-0.7,  albedo 0.3-0.4
                PARA.ensemble.soil.externalWaterFlux(n) = 0.;
                PARA.ensemble.soil.ratioET(n) = 0.;
            case 'shoulder'
                PARA.ensemble.initial_altitude(n) = EBHag - SlopeEB .* (DistanceRC(n)-WidthRoad);
                PARA.ensemble.snow.relative_maxSnow(n) = EBHag - SlopeSnow  .*(DistanceRC(n)-WidthRoad) - PARA.ensemble.initial_altitude(n); % linear deacrease at 0m snow height at the road surface to 40 cm at 5 m distance to the embankment toe
                PARA.ensemble.soil.albedo(n)=0.3;  
                PARA.ensemble.soil.externalWaterFlux(n) = 0.;
                PARA.ensemble.soil.ratioET(n) = 0.;
            case 'toe'
                PARA.ensemble.initial_altitude(n) = 0.;
                PARA.ensemble.snow.relative_maxSnow(n) = EBHag - SlopeSnow  .*(DistanceRC(n)-WidthRoad) - PARA.ensemble.initial_altitude(n); % linear deacrease at 0m snow height at the road surface to 40 cm at 5 m distance to the embankment toe
                PARA.ensemble.soil.albedo(n)=0.2;      % albedo snow-free surface
                PARA.ensemble.soil.externalWaterFlux(n) = 0.002; % external water flux / drainage in [m/day]
                PARA.ensemble.soil.ratioET(n) = 0.5;         
            case 'tundra'
                PARA.ensemble.initial_altitude(n) = 0.;
                PARA.ensemble.snow.relative_maxSnow(n) = rel_maxSnow;  
                PARA.ensemble.soil.albedo(n)=0.2;      % albedo snow-free surface
                PARA.ensemble.soil.externalWaterFlux(n) = 0.002; % external water flux / drainage in [m/day]  temptemp
                PARA.ensemble.soil.ratioET(n) = 0.5;
        end
%% Fuel Storage Tank  
    elseif strcmp(PARA.Exp.Case,'FuelTank')
        %SlopeSnow = (EBHag-rel_maxSnow)/(WidthShoulder+WidthToe);
        if(numlabs>1); SlopeEB = EBHag/sum(PARA.ensemble.TileWidth(3:4)); end% embankment slope  
        switch string(PARA.ensemble.TileType(n))            
            case {'tank_bottom','pile'}
                PARA.ensemble.initial_altitude(n) = EBHag;
                PARA.ensemble.snow.relative_maxSnow(n) = 0.; % max value relative_maxSnow: 0.93m for tile 3 
                PARA.ensemble.soil.albedo(n)=0.3;   % does not apply (no incoming SWR)        
                PARA.ensemble.soil.externalWaterFlux(n) = 0.;
                PARA.ensemble.soil.ratioET(n) = 0.;
                %PARA.soil.evaporationDepth = 0.02; % todotodo  define as ensemble struct member?
            case 'shoulder'    
                PARA.ensemble.initial_altitude(n) = EBHag - SlopeEB .* (DistanceRC(n)-25); % zzz  replace 25 by WidthTank...
               % PARA.ensemble.snow.relative_maxSnow(n) = rel_maxSnow; % todotodo  update...
                PARA.ensemble.snow.relative_maxSnow(3) = 0.5; % zzz
                PARA.ensemble.snow.relative_maxSnow(4) = 1.5; % zzz
                PARA.ensemble.soil.albedo(n)=0.3;      
                PARA.ensemble.soil.externalWaterFlux(n) = -0.1;  % remove snow meltwater  (m/day)
                PARA.ensemble.soil.ratioET(n) = 0.;
                %PARA.soil.evaporationDepth = 0.02;
            case 'toe'    
                PARA.ensemble.initial_altitude(n) = 0.;
                PARA.ensemble.snow.relative_maxSnow(n) = 2.0;
                PARA.ensemble.soil.albedo(n)=0.3;       
                PARA.ensemble.soil.externalWaterFlux(n) = -0.1; % remove snow meltwater
                PARA.ensemble.soil.ratioET(n) = 0.;
                %PARA.soil.evaporationDepth = 0.02;
            case 'foundation_base'
                PARA.ensemble.initial_altitude(n) = 0.;
                PARA.ensemble.snow.relative_maxSnow(n) = 0.5; 
                PARA.ensemble.soil.albedo(n)=0.3;   
                PARA.ensemble.soil.externalWaterFlux(n) = -0.1;
                PARA.ensemble.soil.ratioET(n) = 0.;
                %PARA.soil.evaporationDepth = 0.02;
            case 'disturbed tundra'
                PARA.ensemble.initial_altitude(n) = 0.;
                PARA.ensemble.snow.relative_maxSnow(n) = rel_maxSnow;  
                %PARA.ensemble.snow.relative_maxSnow(n) = 2.0;  
                PARA.ensemble.soil.albedo(n)=0.2;      % albedo snow-free surface
    %           PARA.ensemble.soil.externalWaterFlux(n) = 0.002; % external water flux / drainage in [m/day]  temptemp
                PARA.ensemble.soil.externalWaterFlux(n) = 0.0; % external water flux / drainage in [m/day]  temptemp
                PARA.ensemble.soil.ratioET(n) = 0.5;
            case 'tundra'
                PARA.ensemble.initial_altitude(n) = 0.;
                PARA.ensemble.snow.relative_maxSnow(n) = rel_maxSnow;  
                %PARA.ensemble.snow.relative_maxSnow(n) = 2.0;  
                PARA.ensemble.soil.albedo(n)=0.2;      % albedo snow-free surface
    %           PARA.ensemble.soil.externalWaterFlux(n) = 0.002; % external water flux / drainage in [m/day]  temptemp
                PARA.ensemble.soil.externalWaterFlux(n) = 0.0; % external water flux / drainage in [m/day]  temptemp
                PARA.ensemble.soil.ratioET(n) = 0.5;
        end
%% 1tile (test) case        
    else  
        PARA.ensemble.initial_altitude(n) = 0.;
        PARA.ensemble.snow.relative_maxSnow(n) = rel_maxSnow;  
        PARA.ensemble.soil.albedo(n)=0.2;      % albedo snow-free surface
        PARA.ensemble.soil.externalWaterFlux(n) = 0.0; % external water flux / drainage in [m/day]  temptemp
        PARA.ensemble.soil.ratioET(n) = 0.5;
    end
end
%%        

if(numlabs>1) %ccc
    PARA.IS.TileType = PARA.ensemble.TileType(index);    
    PARA.snow.relative_maxSnow = PARA.ensemble.snow.relative_maxSnow(index);
%    PARA.soil.relative_maxWater= PARA.ensemble.soil.relative_maxWater(index);
    PARA.soil.albedo = PARA.ensemble.soil.albedo(index);
    PARA.soil.externalWaterFlux = PARA.ensemble.soil.externalWaterFlux(index);
    PARA.soil.ratioET = PARA.ensemble.soil.ratioET(index);
end

%A = double( PARA.ensemble.distanceBetweenPoints > 0 ); % adjacency matrix of the network (auxiliary)
A = diag(ones(numlabs-1,1),-1)+diag(ones(numlabs-1,1),1); % adjacency matrix (determines tile interactions)   todotodo  define distances as vector, not as matrix...
PARA.ensemble.distanceBetweenPoints = zeros(numlabs,numlabs) ;
for i=1:numlabs  
    for j=1:numlabs
        if (i-j)==-1
            PARA.ensemble.distanceBetweenPoints(i,j) = TileDistance(i) * A(i,j);
        elseif(i-j)==1
            PARA.ensemble.distanceBetweenPoints(i,j) = TileDistance(i-1) * A(i,j);
        end
    end
end
% % re-define distance of last tile...
% if numlabs==1
%     PARA.ensemble.distanceBetweenPoints = NaN; 
% elseif numlabs==2
%     PARA.ensemble.distanceBetweenPoints = [0 7.5; 7.5 0]; %ccc what is an appropriate distance here?  
%     PARA.ensemble.distanceBetweenPoints = [0 PARA.ensemble.TileWidth(1); PARA.ensemble.TileWidth(1) 0];
% else
%     PARA.ensemble.distanceBetweenPoints(numlabs,numlabs-1) = PARA.ensemble.distanceBetweenPoints(numlabs-1,numlabs-2);
%     PARA.ensemble.distanceBetweenPoints(numlabs-1,numlabs) = PARA.ensemble.distanceBetweenPoints(numlabs-2,numlabs-1);
%     PARA.ensemble.distanceBetweenPoints(numlabs,numlabs) = 0.;
%     %PARA.ensemble.distanceBetweenPoints = mean(PARA.ensemble.TileWidth).* ( diag(ones(numlabs-1,1),-1)+diag(ones(numlabs-1,1),1) ); % equals the distance between mid-points of tiles
% end
 
%PARA.ensemble.weight = [PARA.ensemble.TileWidth(1)/PARA.ensemble.TileWidth(2) 1]; % weight by tile width 
%w1=PARA.ensemble.TileWidth(1)/sum(PARA.ensemble.TileWidth); w2=PARA.ensemble.TileWidth(2)/sum(PARA.ensemble.TileWidth);
%PARA.ensemble.weight = [w1 w2]
PARA.ensemble.weight = PARA.ensemble.TileWidth;
%PARA.ensemble.area = PARA.ensemble.weight ./ sum(PARA.ensemble.weight) .* area_tot ; % in m^2
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  % initial altitudes now defined further above...
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

% parameters related to HEAT exchange   tsvd IS   ccc only needed for lat ex...
switch string(PARA.Exp.Case)
    case 'GravelRoad'
        contact_length = 1000.;
        %PARA.ensemble.thermal_contact_length = 1000. * ones(1,numlabs);  % arbitrary value, describe length in road direction which is unlimited... (cancels out when dT_dt is calculated!)
        PARA.ensemble.thermal_contact_length = contact_length * (diag(1*ones(numlabs-1,1),-1)+diag(1*ones(numlabs-1,1),1)); 
        PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints;
        % thermal distance calculation is now updated in calculateLateralHeatFluxes.m !!!   todotodo make more consistent...
        PARA.ensemble.area = PARA.ensemble.thermal_contact_length .* PARA.ensemble.TileWidth;
    case 'FuelTank'  % ZZZ make consistent with GravelRoad case...use vector in both cases!?
        %PARA.ensemble.thermal_contact_length = 2*pi * PARA.ensemble.TileWidth;  
        PARA.ensemble.thermal_contact_length = 2*pi * cumsum(PARA.ensemble.TileWidth);
        PARA.ensemble.area = [pi*PARA.ensemble.TileWidth(1)^2   pi*( (cumsum(PARA.ensemble.TileWidth(2:numlabs))+PARA.ensemble.TileWidth(1)).^2 - cumsum(PARA.ensemble.TileWidth(1:numlabs-1)).^2)]; % area is defined by plates, besides innermost tile which is a circle
        PARA.ensemble.thermalDistance = TileDistance;
    otherwise % test case   
        PARA.ensemble.thermal_contact_length = 1000.; % dummy values
        PARA.ensemble.area = 10.; 
        PARA.ensemble.thermalDistance = 1000.;
end
       
% parameters related to WATER exchange
%PARA.ensemble.area = [PARA.ensemble.TileWidth(1)*max(PARA.ensemble.thermal_contact_length(:))  PARA.ensemble.TileWidth(2)*max(PARA.ensemble.thermal_contact_length(:))]; % only used for water fluxes...   tsvd IS
%PARA.ensemble.area = max(PARA.ensemble.thermal_contact_length(:)) * [PARA.ensemble.TileWidth(1) PARA.ensemble.TileWidth(2) PARA.ensemble.TileWidth(3) PARA.ensemble.TileWidth(4)]; % only used for water fluxes... cccc  tsvd IS

PARA.ensemble.water_fluxes = zeros( numlabs, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
%tsvd IS PARA.ensemble.water_fluxes = zeros( 1, numlabs ); % total water flux in [m] per sync interval from each worker to worker index
%PARA.ensemble.external_water_flux= zeros( 1, numlabs) ;	%in m/day
% A_scaled = A; A_scaled(1:PARA.IS.number_tiles_IS,1:PARA.IS.number_tiles_IS) = A_scaled(1:PARA.IS.number_tiles_IS,1:PARA.IS.number_tiles_IS)*10. % increase hydraulic conductivity by factor 10 for IS tiles
PARA.ensemble.hydraulic_conductivity = PARA.soil.hydraulic_conductivity * A;
%PARA.ensemble.hydraulic_conductivity = PARA.soil.hydraulic_conductivity * A_scaled;

PARA.ensemble.water_table_altitude = nan(1, numlabs);
%PARA.ensemble.hydraulic_contact_length = PARA.ensemble.thermal_contact_length;
PARA.ensemble.hydraulic_contact_length = 1000 * (diag(1*ones(numlabs-1,1),-1)+diag(1*ones(numlabs-1,1),1));  %temptemp temporal fix!
PARA.ensemble.infiltration_altitude = nan(1, numlabs);
PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints;

%tsvd IS  boundaryCondition={'NoBC','NoBC','NoBC','NoBC'}; 		% set to 'DarcyReservoir' for an external water reservoir
for i=1:numlabs
    boundaryCondition{i}='NoBC';
end
boundaryCondition=transpose(boundaryCondition);
% set to 'DarcyReservoir' for an external water reservoir
%boundaryCondition={'NoBC','NoBC'}; 	% assume that tundra tile can be flooded  ccc
Darcy_elevation=nan(1, numlabs); 	    % Elevation of the Darcy reservoir that can drain or refill the worker it is connected to. NaN for workers without this boundary condition
Darcy_fluxFactor=nan(1, numlabs); 		% Taken as the hydraulic_contact_length*hydraulic_conductivity/hydraulic_distance. Defined for now like this, lets see if we wantto define it differently. NaN for workers without this boundary condition
%tsvd IS   assume that outermost tundra tile is connected to drainage network
if(PARA.soil.drainage==1)
    boundaryCondition{end}='DarcyReservoir';
    Darcy_elevation(end)=0.1; % drainage within the upper 10cm of soil
    Darcy_fluxFactor(end)=PARA.soil.hydraulic_conductivity;
end
PARA.ensemble.boundaryCondition(length(boundaryCondition)).type=boundaryCondition{end};
[PARA.ensemble.boundaryCondition.type]=boundaryCondition{:};
for i=1:numlabs
    if strcmp(boundaryCondition{i},'DarcyReservoir')==1
        PARA.ensemble.boundaryCondition(i).parameters.elevation=Darcy_elevation(i);  
        PARA.ensemble.boundaryCondition(i).parameters.fluxFactor=Darcy_fluxFactor(i); 
    end
end

% parameters related to snow exchange
%tsvd temptemp ZZZ PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight); %ccc weight dependency ok? currently no lat snow!
% if(numlabs>1) %temptemp fix
%     PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, [1 1 1 1 1 1 1]); %ccc weight dependency ok? 
% else
%     PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight);
% end
PARA.ensemble.terrain_index_snow=zeros(1,numlabs);
if(PARA.modules.exchange_snow) %tsvd IS
    PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.altitude, PARA.ensemble.weight); %ccc weight dependency ok? currently no lat snow!
end
%PARA.ensemble.immobile_snow_height = [ 0.1, 0.1, 0.1 ];
PARA.ensemble.immobile_snow_height = 0.1*ones(1,numlabs); % zzz need to set to zero for road surface, or not needed anyhow?
PARA.ensemble.snow_scaling = ones(1, numlabs);  % unclear if needed in ensemble struct

% parameters related to infiltration scheme
% to be specificed by user
PARA.ensemble.rootDepth = 0.2 .* ones(1,numlabs);
%tsvd PARA.ensemble.fieldCapacity = 0.5 .* ones(1, numlabs); % not needed

% location-specific fix parameter values
PARA.location.area = PARA.ensemble.area(index);
PARA.location.initial_altitude = PARA.ensemble.initial_altitude(index);
% location-specific dynamic auxiliary variables
%tsvd PARA.location.area = PARA.ensemble.area(index);
%PARA.location.area = xxx * PARA.ensemble.thermal_contact_length;

PARA.location.altitude = PARA.ensemble.altitude(index);
PARA.location.surface_altitude = PARA.ensemble.surface_altitude(index);
PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(index);
PARA.location.infiltration_altitude = PARA.ensemble.infiltration_altitude(index); 
PARA.location.soil_altitude = PARA.ensemble.soil_altitude(index);

% location-specific dynamic common thresholds
%tsvd IS    use min of soil altitudes for maxWater
%PARA.location.absolute_maxWater_altitude = [max( PARA.ensemble.altitude )+ PARA.soil.relative_maxWater]; % org implementation

%PARA.location.absolute_maxWater_altitude = [min( PARA.ensemble.soil_altitude ) + PARA.soil.relative_maxWater]; % use ensemble.soil_altitude instead of PARA.ensemble.altitude
PARA.location.absolute_maxWater_altitude = PARA.ensemble.altitude(end) + PARA.soil.relative_maxWater; % used fixed value for all tiles, definded w.r.t. last tundra tile

%% soil stratigraphy
zOSL = PARA.Exp.zOSL;  % depth organic surface layer
if strcmp(PARA.Exp.Case,'GravelRoad')
    % column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer
    % extends until the end of the model domain
    % column 2: volumetric water+ice content
    % column 3: volumetric mineral content
    % column 4: volumetric organic content
    % column 5: code for soil type: 1: sand, 2: silt
    % column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs
    BR = 10.0;  % depth bedrock
    switch string(PARA.IS.TileType)
        case 'tundra'
            disp('Gravel Road setting') 
            PARA.soil.layer_properties = ...
              [0.0    0.7     0.05    0.2    4   0.75 ;...   % peat surface layer
               zOSL   0.4     0.55    0.05   2   0.4 ;...    % silty soil (2m thick)
            zOSL+2.0  0.4     0.55    0.05   1   0.4 ;...    % mineral (sandy) soil
               BR     0.3     0.7     0.00   1   0.3 ];      % bedrock                 
            check_layer_properties(PARA) 
        case 'toe'
            if(PARA.modules.xice) % xice at toe (ice wedge / seggregated ice, see SI Raynolds et al. for data)
                PARA.soil.layer_properties = ...
                  [0.0    0.7     0.05    0.2    4   0.75 ;...    % peat surface layer
                   zOSL   0.4     0.55    0.05   2   0.4 ;...    % silty soil - Andersland, Landanyi 1994, table 2.1 silty sand: porosity 23-47%
                   1.0    0.9     0.1     0.0    2   0.4; ...    % xice (50%)  
                   2.0    0.4     0.55    0.05   2   0.4 ;...    % silty soil 
              zOSL+2.0    0.4     0.55    0.05   1   0.4 ;...    % mineral (sandy) soil
                   BR     0.3     0.70    0.00   1   0.3 ];      % bedrock
            else
                PARA.soil.layer_properties = ...
              [0.0    0.7     0.05    0.2    4   0.75 ;...   % peat surface layer
               zOSL   0.4     0.55    0.05   2   0.4 ;...    % silty soil (2m thick)
            zOSL+2.0  0.4     0.55    0.05   1   0.4 ;...    % mineral (sandy) soil
               BR     0.3     0.7     0.00   1   0.3 ];      % bedrock                 
            check_layer_properties(PARA) 
            end
        case 'road' % road center with surface gravel layer 
            PARA.soil.layer_properties = ... 
              [0.0     0.05    0.8      0.0   6    0.2 ;...     % embankment centre - surface layer fine grained(10cm)
               0.1     0.2     0.7      0.0   5    0.3 ;...     % embankment - coarse grained 
               EBT     0.4     0.55     0.05  2    0.4 ;...     % silt (2m)     for porosity values see Andersland, Landanyi 1994, table 2.1 
        EBHag+zOSL+2.0 0.4     0.55     0.05  1    0.4 ;...     % mineral soil
             EBHag+BR  0.3     0.7      0.0   1    0.3 ];       % bedrock   
            check_layer_properties(PARA)
        case 'shoulder' % embankment shoulder
            EBHag_s = PARA.ensemble.initial_altitude(index);       % % embankment height shoulder above ground
            EBT_s =  PARA.ensemble.initial_altitude(index) + EBHbg; % embankment thickness shoulder

            PARA.soil.layer_properties = ...                    % embankment shoulder (no soil cover)
              [0.0       0.2     0.7     0.0    5   0.3 ;...    % embankment (coarse grained)
               EBT_s     0.4     0.55    0.05   2   0.4 ;...    % silt (2m)
        EBHag_s+zOSL+2.0 0.4     0.55    0.05   1   0.4 ;...    % mineral soil
          EBHag_s+BR     0.3     0.7     0.0    1   0.3 ]; ...  % bedrock    
            check_layer_properties(PARA)  
    end
%%  Fuel Storage Tank
elseif strcmp(PARA.Exp.Case,'FuelTank')
    % column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer
    % extends until the end of the model domain
    % column 2: volumetric water+ice content
    % column 3: volumetric mineral content
    % column 4: volumetric organic content
    % column 5: code for soil type: 1: sand, 2: silt
    % column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs
    % 7 concrete
    % 8 steel ???
    BR = 10.0; % depth bedrock
    switch string(PARA.IS.TileType)
        case 'tundra'
            PARA.soil.layer_properties = ...
            [0.0    0.7     0.05    0.2    4   0.75 ;...   % peat surface layer
            zOSL    0.4     0.55    0.05   1   0.4 ;...    % mineral (sandy) soil
            BR      0.3     0.7     0.00   1   0.3 ];      % bedrock    
            check_layer_properties(PARA)        
        case 'disturbed tundra'
            PARA.soil.layer_properties = ...
            [0.0    0.4     0.45    0.05   1   0.5 ;...   % peat surface layer
            zOSL    0.4     0.55    0.05   1   0.4 ;...    % mineral (sandy) soil
            BR      0.3     0.7     0.00   1   0.3 ];      % bedrock    
            check_layer_properties(PARA)        
        case 'pile'
            PARA.soil.layer_properties = ...
            [0.0     0.0     1.0      0.0   8    0.0 ;...    % steel pile 
           EBHag+5.  0.4     0.55    0.05   1    0.4 ;...     % mineral (sandy) soil
           EBHag+BR  0.3     0.7      0.0   1    0.3 ];      % bedrock   
           % check_layer_properties(PARA)
        case 'tank_bottom'    
           PARA.soil.layer_properties = ...
            [0.0     0.0     0.9     0.0   7    0.1 ;...    % concrete layer (50cm)
             0.5     0.05    0.7      0.0   5    0.3 ;...     % gravel   todotodo  define height thickness of concrete and gravel layers, avoid hard-wired...
           EBHag+0.5 0.3     0.55    0.05   1    0.4 ;...     % mineral (sandy)
           EBHag+BR  0.3     0.7      0.0   1    0.3 ];       % bedrock   
           % check_layer_properties(PARA)
        case 'shoulder'   
            EBHag_s = PARA.ensemble.initial_altitude(index);       % % embankment height shoulder above ground
            EBT_s =  PARA.ensemble.initial_altitude(index) + EBHbg; % embankment thickness shoulder
            PARA.soil.layer_properties = ...
            [0.0     0.0     0.9    0.0    7    0.1 ;...    % concrete  todotodo ...update below...
             0.5     0.05     0.7    0.0    5    0.3 ;...     % gravel   todotodo  define height thickness of concrete and gravel layers, avoid hard-wired...
        EBHag_s+0.5  0.3     0.55    0.05   1    0.4 ;...     % mineral (sandy)
        EBHag_s+BR   0.3     0.7      0.0   1    0.3 ];      % bedrock   
           % check_layer_properties(PARA)
        case {'toe','foundation_base'}
             PARA.soil.layer_properties = ...
            [0.0     0.0     0.9    0.0    7    0.1 ;...     % concrete 
             0.5     0.3     0.55    0.05   1    0.4  ;...     % mineral (sandy)
              BR     0.3     0.7      0.0   1    0.3 ];        % bedrock   
           % check_layer_properties(PARA)     
%         case 'foundation_base'    
%             PARA.soil.layer_properties = ...
%             [0.0     0.0     0.9      0.0   7    0.1 ;...    % concrete 
%              1.0     0.4     0.55    0.05   5    0.4 ;...     % gravel   todotodo  define height thickness of concrete and gravel layers, avoid hard-wired...
%              2.0     0.4     0.55    0.05   1    0.4 ;...     % mineral (sandy)
%               BR     0.3     0.7      0.0   1    0.3 ];      % bedrock   
%            % check_layer_properties(PARA)       
    end
%%   1Tile setting 
else
    BR = 10.0; % depth bedrock
    PARA.soil.layer_properties = ...
    [0.0    0.2     0.3    0.1    1   0.6 ;...   % surface layer
     0.1    0.25    0.55   0.05   1   0.4 ;...    % mineral (sandy) soil
     BR     0.3    0.7    0.00   1   0.3 ];      % bedrock    
    check_layer_properties(PARA)    
end


