% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
%
% Developed by: S. Westermann and M. Langer 2015
%
% -------------------------------------------------------------------------

clear
clc
close all
delete(gcp('nocreate')) % useful to restart from a crash

add_modules;  %adds required modules


% Number or realizations
number_of_realizations=2;
if number_of_realizations>1
    parpool(number_of_realizations);
end

% Name, Forcing and saving directory
run_numberi='April12_tryJanreserv2';
forcingname='samoylov_ERA_obs_fitted_1979_2014_spinup.mat';
saveDir = './runs';

diary(['./runs/' run_numberi '_log.txt'])


spmd
    % Name of the run
    run_number=[run_numberi '_real' num2str(labindex)];
    index=labindex;   %number identifying the process; change this to e.g. 1 for single realization (non-parallel) run
    
    %---------------define input parameters------------------------------------
    % here you provide the ground stratigraphy
    % z     w/i     m       o     type porosity
    % default stratigraphy used in publication:
    PARA.soil.layer_properties=[ 0.0   0.60    0.10    0.15    1   0.75;...
        0.15  0.65    0.3     0.05    2   0.4;...
        0.9   0.65    0.3     0.05    1   0.65;...
        9.0   0.30    0.70    0.00    1   0.30     ];
    % soil stratigraphy
    % column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer
    % extends until the end of the model domain
    % column 2: volumetric water+ice content
    % column 3: volumetric mineral content
    % column 4: volumetric organic content
    % column 5: code for soil type: 1: sand, 2: silt
    % column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs
    
    %------ model parameters --------------------------------------------------
    % parameters related to soil
    PARA.soil.albedo=0.2;       % albedo snow-free surface
    PARA.soil.epsilon=0.97;     % emissvity snow-free surface
    PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
    PARA.soil.rs=50;            % surface resistance against evapotransiration [m^-1] snow-free surface
    PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
    PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]
    
    % parameters related to hydrology scheme
    PARA.soil.fieldCapacity=0.5;    %water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
    PARA.soil.evaporationDepth=0.10; %depth to which evaporation occurs - place on grid cell boundaries
    PARA.soil.rootDepth=0.2;        %depth affected by transpiration - place on grid cell boundaries
    %    PARA.soil.wiltingPoint=0.2;     %point at which transpiration shuts off
    %    PARA.soil.residualWC=0.05;      %water always remaining in the soil, not accessible to evaporation
    PARA.soil.ratioET=0.5;          % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
    PARA.soil.externalWaterFlux=0;  %external water flux / drainage in [m/day]
    PARA.soil.convectiveDomain=[];       % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
    PARA.soil.mobileWaterDomain=[0 10.0];      % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
    PARA.soil.relative_maxWater=0;              % depth at which a water table will form [m] - above excess water is removed, below it pools up
    PARA.soil.hydraulic_conductivity = 1e-5;
    PARA.soil.infiltration_limit_depth=1.25;     % depth [m] from the surface at wich the infiltration bucket scheme is stopped if no permafrost.
    PARA = loadSoilTypes( PARA );
    
    % parameters related to snow
    PARA.snow.max_albedo=0.85;      % albedo of fresh snow
    PARA.snow.min_albedo=0.5;       % albedo of old snow
    PARA.snow.epsilon=0.99;         % surface emissivity snow
    PARA.snow.z0=5e-4;              % roughness length surface [m]
    PARA.snow.rs=0.0;               % surface resistance -> should be 0 for snow
    PARA.snow.rho_snow=200.0;       % density in [kg/m3]
    PARA.snow.tau_1=86400.0;        % time constants of snow albedo change (according to ECMWF reanalysis) [sec]
    PARA.snow.tau_a=0.008;          % [per day]
    PARA.snow.tau_f=0.24;           % [per day]
    PARA.snow.relative_maxSnow= 1.0; 	% maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.snow.extinction=25.0;      % light extinction coefficient of snow
    
    % parameters related to water body on top of soil domain
    PARA.water.albedo=0.07;     % albedo water (parameterization after Wayne and Burt (1954) in surfaceCondition.m)
    PARA.water.epsilon=0.99;    % surface emissivity water
    PARA.water.rs=0.0;            % surface resistance -> should be 0 for water
    PARA.water.z0=5e-4;              % roughness length surface [m] % JAN: value for summer / vegetation
    
    PARA.ice.albedo =0.20;      % albedo ice / Lei et al. (2011) shows a range of 0.1 to 0.35
    PARA.ice.epsilon=0.98;      % surface emissivity snow
    PARA.ice.rs=0.0;              % surface resistance -> should be 0 for ice
    PARA.ice.z0=5e-4;              % roughness length surface [m] % JAN: value for snow
    PARA.ice.extinction=4.5;    % [m^-1] light extinction coefficient of ice / Lei et al. (2011) shows a range of 1 to 5 m^-1
    
    PARA.technical.z=2.0;                       % height of input air temperature above ground in [m] - assumed constant even when snow depth increases
    PARA.technical.SWEperCell=0.005;            % SWE per grid cell in [m] - determines size of snow grid cells
    PARA.technical.maxSWE=0.4;                  % in [m] SWE
    PARA.technical.arraySizeT=5002;             % number of values in the look-up tables for conductivity and capacity
    PARA.technical.starttime=datenum(1979, 7, 1);       % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=datenum(1979, 7, 10);         % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;   % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;   % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;            % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3 ./ 24.0 ;          % output time step in [days] - here three hours
    PARA.technical.syncTimeStep = 12 ./ 24.0 ;          % output time step in [days] - here three hours
    PARA.technical.saveDate='01.01.';           % date of year when output file is written - no effect if "saveInterval" is empty
    PARA.technical.saveInterval=1;             % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;          % default size of a newly added water cell when water ponds below water table [m]
    
    %default grid used for publications and testing of water balance:
    PARA.technical.subsurfaceGrid = [0:0.02:4, 4.1:0.1:10, 10.2:0.2:20, 21:1:30, 35:5:50, 60:10:100, 200:100:1000]'; % the subsurface K-grid in [m]
    %PARA.technical.subsurfaceGrid = [[0:0.02:10], [10.1:0.1:20], [20.2:0.2:30], [31:1:40], [45:5:60], [70:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
    
    PARA.location.area=1.0;
    PARA.location.initial_altitude=20.0;
    % dynamic auxiliary varaibles
    PARA.location.altitude = PARA.location.initial_altitude; 	% used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given; excluding snow domain
    PARA.location.surface_altitude=PARA.location.initial_altitude;		% this is dynamic and refers to the surface including snow
    PARA.location.soil_altitude = PARA.location.initial_altitude;
    PARA.location.infiltration_altitude = nan; % defined at runtime
    PARA.location.water_table_altitude = nan; % defined at runtime
    PARA.soil.infiltration_limit_altitude=PARA.location.initial_altitude-PARA.soil.infiltration_limit_depth; % LEO: for single run, but updated for multiple runs as the min of all workers
    PARA.location.bottomBucketSoilcTIndex = nan; % defined at runtime
    % common thresholds
    PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
    if isempty( PARA.snow.relative_maxSnow )
        PARA.location.absolute_maxSnow_altitude = [];
    else
        PARA.location.absolute_maxSnow_altitude = PARA.location.altitude + PARA.snow.relative_maxSnow;
    end
    
    %initial temperature profile -> first column depth [m] -> second column temperature [degree C]
    PARA.Tinitial = [  -2     5   ;...
        0     0   ;...
        2    -2   ;...
        5    -7   ;...
        10    -9  ;...
        25    -9   ;...
        100    -8   ;...
        1100    10.2   ];      % the geothermal gradient for Qgeo=0.05W/m^2 and K=2.746W/Km is about 18.2 K/km
    
    
    PARA = loadConstants( PARA );
    
    %FORCING data mat-file
    PARA.forcing.filename=forcingname;  %must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
    PARA.forcing.rain_fraction=1;
    PARA.forcing.snow_fraction=1;
    
    % switches for modules
    PARA.modules.infiltration=1;   % true if infiltration into unfrozen ground occurs
    PARA.modules.xice=1;           % true if thaw subsicdence is enabled
    PARA.modules.lateral=1;		   % true if adjacent realizations are run (this does not require actual lateral fluxes)
    
    if PARA.modules.lateral
        % switches for lateral processes
        PARA.modules.exchange_heat = 1;
        PARA.modules.exchange_water = 1;
        PARA.modules.exchange_snow = 1;
        
        %---------overwrites variables for each realization--------------------
        % this function must define everything that is realization-specific or dependent of all realizations
        PARA = get_parallel_variables( PARA );
    end
    
    % ------make output directory (name depends on parameters) ----------------    
    mkdir([ saveDir '/' run_number]);
    
    %--------------------------------------------------------------------------
    %-----------do not modify from here onwards--------------------------------
    %--------------------------------------------------------------------------
    [FORCING, success]=load_forcing_from_file(PARA); % load FORCING mat-file
    
    if ~success
        warning('A problem with the Forcing occured.');
    end
    
    PARA = initializeParameters(PARA, FORCING); %set start time, etc.
    
    %----------------create and initialize the grids --------------------------
    GRID=makeGrids(PARA);  %create all grids
    GRID=createStratigraphy(PARA,GRID);   %interpolate input stratigraphy to the soil grid
    
    %----- initializie excess ground ice --------------------------------------
    [GRID,PARA] = initializeExcessIce2(GRID,PARA);
    
    %----- initializie soil thermal properties --------------------------------
    GRID = initializeSoilThermalProperties(GRID, PARA);
    
    %------ initializie snow properties----------------------------------------
    GRID = initializeSnow(GRID);
    
    %---- initialize the surface energy balance struct ------------------------
    SEB = initializeSEB();
    
    %---- initialize the water body module ------------------------------------
    GRID = initializeLAKE(GRID);
    
    %---- initialize temperature profile --------------------------------------
    T = inititializeTemperatureProfile_simple(GRID, PARA, FORCING);
    
    %---- modification for infiltration
    wc=GRID.soil.cT_water;
    GRID.soil.E_lb = find(PARA.soil.evaporationDepth==GRID.soil.soilGrid(:,1))-1;
    GRID.soil.T_lb= find(PARA.soil.rootDepth==GRID.soil.soilGrid(:,1))-1;
    GRID.soil.water2pool=0; % Leo : cannot find a good initialize function to put it in
    
    %---- preallocate temporary arrays for capacity and conductivity-----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA); % this is basically the same as "getThermalProperties" during integration, but without interpolation to K grid
    
    %---- energy and water balance initialization -----------------------------
    BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);
    
    %---- temporary arrays for storage of lateral fluxes --> these could go into a LATERAL struct or TEMPORARY?
    %water_fluxes = zeros( 1, numlabs );                         % total water flux in [m/s] from each worker to worker index
    %snow_fluxes = zeros( 1, numlabs );                          % total snow flux in [m SWE] per output interval from each worker to worker index
    %heat_fluxes = zeros( 1, numlabs );                           % total heat flux in [J] per output interval of all workers to worker index
    %dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s]
    
    %__________________________________________________________________________
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    
    fprintf('initialization successful\n');
    iSaveSettings(  [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_settings.mat'] , FORCING, PARA, GRID)
    
    
    %% ________________________________________________________________________
    % Time Integration Routine                                                I
    %                                                                         I
    %_________________________________________________________________________I
    
    while t<PARA.technical.endtime
        
        %------ interpolate forcing data to time t ----------------------------
        FORCING = interpolateForcingData(t, FORCING);
        
        %------determine the thermal properties of the model domains ----------
        [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = getThermalPropertiesInfiltration(T, wc, c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid, GRID, PARA);
        
        %------- water and energy balance calculations ------------------------
        BALANCE = updateBALANCE(T, wc, c_cTgrid, lwc_cTgrid, BALANCE, GRID, PARA);
        
        %------ surface energy balance module ---------------------------------
        %set surface conditions (albedo, roughness length, etc.)
        [PARA, GRID] = surfaceCondition(GRID, PARA, T);
        %calculate the surface energy balance
        [SEB, dwc_dt] = surfaceEnergyBalanceInfiltration(T, wc, FORCING, GRID, PARA, SEB);
        
        %------ soil module  --------------------------------------------------
        %calculate heat conduction
        SEB = heatConduction(T, k_Kgrid, GRID, PARA, SEB);
        
        %------ sum up heat fluxes --------------------------------------------
        SEB.dE_dt = SEB.dE_dt_cond + SEB.dE_dt_SEB;
        
        %------ determine optimal timestep ------------------------------------
        % account for min and max timesteps specified, max. energy change per grid cell and the CFT stability criterion.
        % energy change due to advection of heat through water fluxes is still excluded.
        % timestep in [days]
        timestep = min( [ max( [ min( [ 0.5 * nanmin( GRID.general.K_delta.^2 .* c_cTgrid ./ k_cTgrid ./ (GRID.soil.cT_domain + GRID.snow.cT_domain ) ) ./ (24.*3600), ...
            PARA.technical.targetDeltaE .* nanmin( abs(GRID.general.K_delta ./ SEB.dE_dt ) ) ./ (24.*3600), ...
            PARA.technical.maxTimestep ] ), ...
            PARA.technical.minTimestep ] ), ...
            TEMPORARY.outputTime-t ] );
        if PARA.modules.lateral
            timestep = min( [timestep, TEMPORARY.syncTime-t] );
        end
        
        % give a warning when timestep required by CFT criterion is below the minimum timestep specified
        if timestep > 0.5 * min( GRID.general.K_delta.^2 .* c_cTgrid ./ k_cTgrid ./ (GRID.soil.cT_domain + GRID.snow.cT_domain) ) / (24.*3600)
            warning( 'numerical stability not guaranteed' );
        end
        
        %------ update T array ------------------------------------------------
        % account for vertical heat fluxes from ground heat flux and heat conduction
        T = T + SEB.dE_dt./c_cTgrid./GRID.general.K_delta.*timestep.*24.*3600;
        % account for lateral heat fluxes --> this is now done directly at sync time
        %T = T + dE_dt_lateral./c_cTgrid.*timestep.*24.*3600; % no division by K_delta necessary as dE_dt_lateral in [ J / m^3 / s ]
        % set grid cells in air to air temperature
        T(GRID.air.cT_domain)=FORCING.i.Tair;
        
        %------- water body module --------------------------------------------
        T = mixingWaterBody(T, GRID);
        
        %------- snow cover module --------------------------------------------
        [T, GRID, PARA, SEB, BALANCE] = CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_cTgrid, timestep, BALANCE);
        [GRID, T, BALANCE] = updateGRID_snow(T, GRID, PARA, BALANCE);
        
        %------- infiltration module-------------------------------------------
        if PARA.modules.infiltration
            [wc, GRID, BALANCE] = CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING, BALANCE);
        end
        
        %------- excess ice module --------------------------------------------
        if PARA.modules.xice && ~PARA.modules.infiltration
            warning( 'energy and water balances are not correct for this combination of modules');
            [GRID, PARA] = excessGroundIce(T, GRID, PARA);
            % assure wc has correct length
            wc = wc( end-sum(GRID.soil.cT_domain)+1 : end );
        elseif PARA.modules.xice && PARA.modules.infiltration
            [GRID, PARA, wc, meltwaterGroundIce] = excessGroundIceInfiltration(T, wc, GRID, PARA);
            GRID = updateGRID_excessiceInfiltration2(meltwaterGroundIce, GRID);
        end
        
        %------- update Lstar for next time step ------------------------------
        SEB = L_star(FORCING, PARA, SEB);
        
        %------- update auxiliary state variables
        PARA.location.altitude = getAltitude( PARA, GRID );
        PARA.location.surface_altitude = getSurfaceAltitude( PARA, GRID );
        PARA.location.soil_altitude = getSoilAltitude( PARA, GRID );
        [PARA.location.infiltration_altitude, PARA.location.bottomBucketSoilcTIndex] = getInfiltrationAltitude( PARA, GRID, T);
        PARA.location.water_table_altitude = getWaterTableAltitudeFC(T, wc, GRID, PARA);
        PARA.soil.infiltration_limit_altitude = PARA.location.soil_altitude - PARA.soil.infiltration_limit_depth;
        
        %------- update threshold variables if no lateral exchange processes occur, otherwise updated at sync time
        if ~PARA.modules.lateral
            PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
            if isempty( PARA.snow.relative_maxSnow ) % Leo : changed from PARA.snow.maxSnow
                PARA.location.absolute_maxSnow_altitude = [];
            else
                PARA.location.absolute_maxSnow_altitude =  PARA.location.altitude + PARA.snow.relative_maxSnow; % Leo : changed ensemble to location otherwise it craches.
            end
        end
        
        %------- water balance calculations -----------------------------------
        % rainfall
        BALANCE.water.dp_rain = BALANCE.water.dp_rain + FORCING.i.rainfall.*timestep;   %sum up rainfall in [mm] per output interval
        % snowfall
        BALANCE.water.dp_snow = BALANCE.water.dp_snow + FORCING.i.snowfall.*timestep;   %sum up snowfall in [mm] SWE per output interval
        
        %------- lateral exchange module --------------------------------------
        % all functions called in this block should go into
        % /modules/cryoGridLateral
        % calling PARA.ensemble is only allowed here
        if PARA.modules.lateral
            if t==TEMPORARY.syncTime %communication between workers
                fprintf('\n\t\t\tCryoGridLateral: sync - start (Worker %1.0f)\n', labindex);
                labBarrier(); %common start
                PARA = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA) ;
                
                % heat exchange module
                if PARA.modules.exchange_heat
                    labBarrier();
                    % check preconditions
                    precondition_heatExchange = true; %no specific conditions so far
                    if precondition_heatExchange
                        %WRAPPER
                        fprintf('\t\t\tsync - exchanging heat\n');
                        % calculate lateral heat fluxes
                        PACKAGE_heatExchange.T = T;
                        PACKAGE_heatExchange.cT_grid = GRID.general.cT_grid;
                        PACKAGE_heatExchange.k_cTgrid = k_cTgrid;
                        for j=1:number_of_realizations
                            if j~=index
                                labSend( PACKAGE_heatExchange, j, 1);
                            end
                        end
                        % provide temporary array for lateral heat fluxes
                        dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s]
                        for j=1:number_of_realizations
                            if j~=index
                                PACKAGE_heatExchange_j = labReceive(j, 1);
                                [dE_dt_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/m^3/s ]
                                TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + dE_dt_lateral_j .* PARA.technical.syncTimeStep.*24.*3600; % in [ J/m^3 ]
                                TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  dE_dt_lateral_j .* PARA.technical.syncTimeStep.*24.*3600 .* GRID.general.K_delta ); % depth intergrated in [ J/m^2 ]
                                dE_dt_lateral = dE_dt_lateral + dE_dt_lateral_j;    % sum up contributions from all realizations in [ J/m^3/s ]
                            end
                        end
                        % apply lateral heat fluxes for entire sync interval directly
                        T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600; % no division by K_delta necessary as dE_dt_lateral in [ J / m^3 / s ]
                    end
                end
                
                % water exchange module
                if PARA.modules.exchange_water
                    labBarrier();
                    % check preconditions
                    precondition_waterExchange = checkPreconditionWaterExchange( T, GRID );
                    if precondition_waterExchange
                        % WRAPPER
                        fprintf('\t\t\tsync - exchanging water\n');
                        % calculate lateral water fluxes
                        water_fluxes= zeros(numlabs,numlabs); % in m of height change
                        PACKAGE_waterExchange.water_table_altitude = PARA.ensemble.water_table_altitude(index);     % JAN: I think we do not need to exchange water_table and infiltration altitude because this is done already at the beginning of the lateral scheme
                        PACKAGE_waterExchange.infiltration_altitude = PARA.ensemble.infiltration_altitude(index);
                        PACKAGE_waterExchange.infiltration_condition = T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub);
                        
                        for j=1:number_of_realizations
                            if j~=index
                                labSend( PACKAGE_waterExchange, j, 2);
                            end
                        end
                        for j=1:number_of_realizations
                            if j~=index
                                PACKAGE_waterExchange_j = labReceive(j, 2);
                                water_fluxes = calculateLateralWaterDarcyFluxes( T, PACKAGE_waterExchange_j, GRID, PARA, j, water_fluxes);  % matrix containing all fluxes in [m/s] scaled to row index
                            end
                        end
                        
                        % Calculate possible boundary fluxes
                        [ boundary_water_flux ] = calculateLateralWaterBoundaryFluxes(PARA, GRID, T);
                        
                        % Check for water availability and set real water fluxes
                        [ water_fluxes_worker, boundary_water_flux ] = calculateLateralWaterAvailable( PARA,GRID, wc, water_fluxes,boundary_water_flux );
                        water_fluxes_gather=zeros(numlabs,numlabs,numlabs);
                        water_fluxes_gather(:,:,labindex)=water_fluxes_worker;
                        
                        % Send real fluxes all around
                        for j=1:number_of_realizations
                            if j~=index
                                labSend( water_fluxes_worker, j, 5);
                            end
                        end
                        
                        for j=1:number_of_realizations
                            if j~=index
                                water_fluxes_worker_j = labReceive(j, 5);
                                water_fluxes_gather(:,:,j)=water_fluxes_worker_j;
                            end
                        end
                        
                        water_fluxes=nansum(water_fluxes_gather,3);
                        waterflux=nansum(water_fluxes(labindex,:))+boundary_water_flux;
                        
                        % apply lateral water flux directly (as bulk subsurface flux)
                        [wc, excess_water, lacking_water] = bucketScheme(T, wc, zeros( size(wc) ), GRID, PARA, waterflux);
                        try
                            assert( lacking_water < 1e-9, 'CryoGrid3 - lateral exchange - lacking water>0');    % there should be no lacking water as this was checked for
                        catch
                            fprintf(' Lacking water = %f m\n', lacking_water );
                        end
                        
                        % Store and display
                        BALANCE.water.dr_water_fluxes_out=BALANCE.water.dr_water_fluxes_out+water_fluxes./(PARA.technical.syncTimeStep*24*3600); % here we decide in which units we want it. Now m/s
                        BALANCE.water.dr_lateralWater = BALANCE.water.dr_lateralWater + (waterflux-excess_water)*1000; % Excess water is removed so that we only keep the net water modification implied by the lateral fluxes
                        fprintf('\t\t\tNet wc change :\t%3.2e m\n',waterflux-excess_water)
                        if excess_water>1e-9
                            GRID.soil.water2pool= GRID.soil.water2pool + excess_water;
                            fprintf('\t\t\tExcess water :\t%3.2e m\n',excess_water)
                            BALANCE.water.dr_lateralExcess=BALANCE.water.dr_lateralExcess + excess_water*1000; % Added by Leo to have the lateral fluxes in BALANCE
                        end
                        if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1;
                            BALANCE.water.dr_DarcyReservoir = BALANCE.water.dr_DarcyReservoir + boundary_water_flux*1000;
                            fprintf('\t\t\tBoundary contribution :\t%3.2e m\n',boundary_water_flux)
                        end
                    end
                    
                end
                
                % snow exchange module
                if PARA.modules.exchange_snow
                    labBarrier();
                    % check preconditions
                    precondition_snowExchange = checkPreconditionSnowExchange( GRID, PARA );
                    if precondition_snowExchange
                        fprintf('\t\t\tsync - exchanging snow\n');
                        % calculate terrain index with updated surface_altitudes
                        PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.surface_altitude, PARA.ensemble.weight);
                        % calculate mobile snow
                        % WRAPPER
                        mobile_snow = zeros( 1, number_of_realizations );
                        my_mobile_snow = 0;
                        meltingConditions_index = sum(GRID.snow.Snow_w)>0 || sum(T(GRID.snow.cT_domain)>0)>0;    %snow is assumed to be immobile under melting conditions CHECK WITH SNOW MODULE
                        if ~isempty(GRID.snow.cT_domain_ub) && ~meltingConditions_index	% current realization has snow cover and no melting conditions
                            i=0;
                            while (abs( GRID.general.K_grid(GRID.snow.cT_domain_ub+i)-GRID.general.K_grid(GRID.snow.cT_domain_lb+1) )... 	% snow only mobile above realization-specific threshold
                                    > PARA.ensemble.immobile_snow_height(index)+GRID.general.K_delta(GRID.snow.cT_domain_ub)) ...           % if upper cell is drifted away, immobile snow height remains
                                    && (PARA.ensemble.initial_altitude(index)-GRID.general.K_grid(GRID.snow.cT_domain_ub+i) ...					% snow only mobile above lowermost surface altitude + snowCellSize (to prevent oscillations)
                                    - (min( PARA.ensemble.surface_altitude )+GRID.snow.snowCellSize) > 1e-6 )
                                my_mobile_snow = my_mobile_snow + (GRID.snow.Snow_i(GRID.snow.cT_domain_ub+i)); %only "ice" mobile
                                i=i+1;
                            end
                        end
                        mobile_snow(index) = my_mobile_snow;
                        % exchange mobile snow amounts
                        for j=1:number_of_realizations
                            if j~=index
                                % send mobile snow amount in [m SWE]
                                labSend( mobile_snow(index), j, 4 );
                            end
                        end
                        for j=1:number_of_realizations
                            if j~=index
                                % receive mobile snow amount [m SWE]
                                mobile_snow(j) = labReceive(j, 4);
                            end
                        end
                        % calculate lateral snow fluxes
                        my_snow_change = calculateLateralSnowFluxes2( mobile_snow, PARA );
                        % apply lateral snow fluxes directly
                        if my_snow_change ~= 0
                            [T, GRID] = applyLateralSnowFluxes( T, PARA, GRID, FORCING, my_snow_change );
                            [GRID, T, BALANCE] = updateGRID_snow(T, GRID, PARA, BALANCE);
                            BALANCE.water.dr_lateralSnow = BALANCE.water.dr_lateralSnow + my_snow_change * 1000;
                        end
                        
                        TEMPORARY.snow_flux_lateral = TEMPORARY.snow_flux_lateral + my_snow_change;
                        fprintf('Snow flux to worker %d = %f mm SWE \n', [index, my_snow_change*1000] );
                        
                    end
                end
                
                labBarrier();
                PARA = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA) ;
                
                TEMPORARY.syncTime=round((TEMPORARY.syncTime + PARA.technical.syncTimeStep)./PARA.technical.syncTimeStep).*PARA.technical.syncTimeStep;
                fprintf('sync - done\n');
            end
        end
        
        
        %------- next time step -----------------------------------------------
        t=t+timestep;
        %---------- sum up + OUTPUT -------------------------------------------
        [TEMPORARY, OUT, BALANCE] = sum_up_output_store(t, T, wc, lwc_cTgrid(GRID.soil.cT_domain), timestep, TEMPORARY, BALANCE, PARA, GRID, SEB, OUT, saveDir, run_number);
        
    end
    
    % save final state and output at t=endtime
    iSaveOUT( [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_output' datestr(t,'yyyy') '.mat'], OUT)
    iSaveState( [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_finalState' datestr(t,'yyyy') '.mat'], T, wc, t, SEB, PARA, GRID)
end

if number_of_realizations>1
    delete(gcp('nocreate'))
end

fprintf('Done with %s\n', run_numberi);
diary off