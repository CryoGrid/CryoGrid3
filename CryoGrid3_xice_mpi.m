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

% dbstop if error;

number_of_realizations=2;

if number_of_realizations>1
    parpool(number_of_realizations);
end

% Name and Forcing
run_number='debug06March_5yrs';
forcingname='Suossjavri_WRF_Norstore_adapted5yr.mat';

spmd
    % Name of the run
    run_number=[run_number '_real' num2str(labindex)];
     
    index=labindex;   %number identifying the process; change this to e.g. 1 for single realization (non-parallel) run
    
    %---------------define input parameters------------------------------------
    
    % default stratigraphy used in publication:
    PARA.soil.layer_properties=[ 0.0   0.55    0.05    0.15   1   0.80;...
                                 0.5   0.80    0.05    0.15   2   0.80;...
                                 3.0   0.50    0.50    0.0    1   0.50;...
                                10.0   0.03    0.97    0.0    1   0.03 ];

    % soil stratigraphy
    % column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer extends until the end of the model domain
    % column 2: volumetric water+ice content
    % column 3: volumetric mineral content
    % column 4: volumetric organic content
    % column 5: code for soil type: 1: sand, 2: silt
    % column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs
    
    %------ model parameters --------------------------------------------------
    % parameters related to soil
    PARA.soil.albedo=0.2;       % albedo snow-free surface
    PARA.soil.albedoPond=0.07;  % albedo of water, used when the uppermost grod cell is 100% water due to modeled thermokarst development
    PARA.soil.epsilon=0.97;     % emissvity snow-free surface
    PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
    PARA.soil.rs=0;             % surface resistance against evapotransiration [m^-1] snow-free surface
    PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
    PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]
    
    % parameters related to hydrology scheme
    PARA.soil.fieldCapacity=0.55;       % water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
    PARA.soil.evaporationDepth=0.05;    % depth to which evaporation occurs - place on grid cell boundaries
    PARA.soil.rootDepth=0.2;            % depth affected by transpiration - place on grid cell boundaries
    PARA.soil.wiltingPoint=0.2;         % point at which transpiration shuts off
    PARA.soil.residualWC=0.05;          % water always remaining in the soil, not accessible to evaporation
    PARA.soil.ratioET=0.5;              % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
    
    
    PARA.soil.externalWaterFlux=0;       % external water flux / drainage in [m/day] % LEO : Parameters to check
    PARA.soil.convectiveDomain=[];       % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
    PARA.soil.mobileWaterDomain=[];      % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile % LEO : Parameters to check
    PARA.soil.relative_maxWater=0;              % depth at which a water table will form [m] - above excess water is removed, below it pools up % LEO : Parameters to check
    PARA.soil.hydraulic_conductivity = 1e-5;
    PARA.soil.infiltration_limit=1.25;     % depth [m] from the surface at wich the infiltration bucket scheme is stopped if no permafrost.
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
    PARA.snow.relative_maxSnow= 0.13; 	% maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.snow.extinction=30.0;      % light extinction coefficient of snow
    
    % parameters related to water body on top of soil domain
    PARA.water.albedo=0.05;     % albedo water (parameterization after Wayne and Burt (1954) in surfaceCondition.m)
    PARA.water.epsilon=0.99;    % surface emissivity water
    PARA.water.rs=0.0;            % surface resistance -> should be 0 for water
    PARA.water.z0=1e-3;              % roughness length surface [m] % JAN: value for summer / vegetation
    
    PARA.ice.albedo =0.20;      % albedo ice / Lei et al. (2011) shows a range of 0.1 to 0.35
    PARA.ice.epsilon=0.98;      % surface emissivity snow
    PARA.ice.rs=0.0;              % surface resistance -> should be 0 for ice
    PARA.ice.z0=5e-4;              % roughness length surface [m] % JAN: value for snow
    PARA.ice.extinction=4.5;    % [m^-1] light extinction coefficient of ice / Lei et al. (2011) shows a range of 1 to 5 m^-1
    
    PARA.technical.z=2.0;                       % height of input air temperature above ground in [m] - assumed constant even when snow depth increases
    PARA.technical.SWEperCell=0.005;            % SWE per grid cell in [m] - determines size of snow grid cells
    PARA.technical.maxSWE=0.4;                  % in [m] SWE
    PARA.technical.arraySizeT=5002;             % number of values in the look-up tables for conductivity and capacity
    PARA.technical.starttime=[]; % datenum(1979, 7, 1);       % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=[];         % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;   % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;   % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;            % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3 / 24;          % output time step in [days] - here three hours
    PARA.technical.syncTimeStep = 12 / 24;          % output time step in [days] - here three hours
    PARA.technical.saveDate='01.08.';           % date of year when output file is written - no effect if "saveInterval" is empty
    PARA.technical.saveInterval=1;              % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;          % default size of a newly added water cell when water ponds below water table [m]
    
    %default grid used for publications and testing of water balance:
    PARA.technical.subsurfaceGrid = [0:0.05:2 2.1:0.5:8 8.2:0.2:20 21:1:30 35:5:50 60:10:100 200:100:1000]'; % the subsurface K-grid in [m]
    %PARA.technical.subsurfaceGrid = [[0:0.02:10], [10.1:0.1:20], [20.2:0.2:30], [31:1:40], [45:5:60], [70:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
    
    PARA.location.area=1.0;
    PARA.location.initial_altitude=300;
    % JAN: the following quantities are dynamic and should hence be moved to another struct, e.g. "STATE"
    PARA.location.altitude = PARA.location.initial_altitude; 	% used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given; excluding snow domain
    PARA.location.surface_altitude=PARA.location.initial_altitude;		% this is dynamic and refers to the surface including snow
    PARA.location.active_layer_depth_altitude = nan; % defined at runtime
    PARA.location.water_table_altitude = nan; % defined at runtime
    PARA.soil.alt_infiltration_limit=PARA.location.initial_altitude-PARA.soil.infiltration_limit; % LEO: for single run, but updated for multiple runs as the min of all workers
    PARA.location.bottomBucketSoilcTIndex = 1;
    % thresholds
    PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
    if isempty( PARA.snow.relative_maxSnow )
        PARA.location.absolute_maxSnow_altitude = [];
    else
        PARA.location.absolute_maxSnow_altitude =  PARA.location.altitude + PARA.snow.relative_maxSnow ;
    end
    
    %initial temperature profile -> first column depth [m] -> second column temperature [degree C]
    PARA.Tinitial = [ -5    10   ;...
                       0     5   ;...
                       0.1   2   ;...
                       0.5   0.5 ;...
                       1     0   ;...
                       2    -0.2 ;...
                      10     0   ;...
                      30     2   ;...
                     500     4   ;...
                    5000    10   ];
    
    PARA = loadConstants( PARA );
    
    %FORCING data mat-file
    PARA.forcing.filename=forcingname;  %must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
    PARA.forcing.rain_fraction=1;
    PARA.forcing.snow_fraction=1;
    
    % switches for modules
    PARA.modules.infiltration=1;   % true if infiltration into unfrozen ground occurs
    PARA.modules.xice=0;           % true if thaw subsicdence is enabled
    PARA.modules.lateral=1;		   % true if adjacent realizations are run (this does not require actual lateral fluxes)
    
    if PARA.modules.lateral
        % switches for lateral processes
        PARA.modules.exchange_heat = 0;
        PARA.modules.exchange_water = 1;
        PARA.modules.exchange_snow = 0;
        
        %---------overwrites variables for each realization--------------------
        % this function must define everything that is realization-specific or dependent of all realizations
        PARA = get_parallel_variables( PARA );
    end
    
    % ------make output directory (name depends on parameters) ----------------
   
    %     run_number= sprintf( 'testrunMPI_POOL_xH%d_xW%d_xS%d_infil%d_xice%d_rF%d_sF%d_realization%d' , ...
    %         [ PARA.modules.exchange_heat, PARA.modules.exchange_water, PARA.modules.exchange_snow, ...
    %         PARA.modules.infiltration, PARA.modules.xice, ...
    %         PARA.forcing.rain_fraction, PARA.forcing.snow_fraction ,  ...
    %         index ] );
    
    mkdir(['./runs/' run_number])
    
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
    
    %---- preallocate temporary arrays for capacity and conductivity-----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA); % this is basically the same as "getThermalProperties" during integration, but without interpolation to K grid
    
    %---- energy and water balance initialization -----------------------------
    BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);
    
    %---- temporary arrays for storage of lateral fluxes --> these could go into a LATERAL struct or TEMPORARY?
    water_fluxes = zeros( 1, numlabs );                         % total water flux in [m/s] from each worker to worker index
    snow_fluxes = zeros( 1, numlabs );                          % total snow flux in [m SWE] per output interval from each worker to worker index
    heat_fluxes = zeros( 1, numlabs );                           % total heat flux in [J] per output interval of all workers to worker index
    dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);    % cell-wise heat flux in [W/m^3]? of all workers to worker index
    
    %__________________________________________________________________________
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    
    disp('initialization successful');
    iSaveSettings( [ './runs/' run_number '/' run_number '_settings.mat'] , FORCING, PARA, GRID)
    
    
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
            TEMPORARY.syncTime-t,...
            TEMPORARY.outputTime-t ] );
        
        
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
            %lateral_water_flux = nansum( water_fluxes );  % lateral fluxes
            %to/from other workers in [m/s] --> this is now done directly
            %at sync time
            [wc, GRID, BALANCE] = CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING, BALANCE, 0);
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
        [bottomBucket_altitude, bottomBucketSoilcTIndex]= getActiveLayerDepthAltitude( PARA, GRID, T);
        PARA.location.active_layer_depth_altitude = bottomBucket_altitude;
        PARA.location.bottomBucketSoilcTIndex = bottomBucketSoilcTIndex;
        PARA.location.water_table_altitude = getWaterTableAltitudeFC(T, wc, GRID, PARA);
        
        %------- update threshold variables if no lateral exchange processes occur, otherwise updated at sync time
        if ~PARA.modules.lateral
            PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
            if isempty( PARA.snow.maxSnow )
                PARA.location.absolute_maxSnow_altitude = [];
            else
                PARA.location.absolute_maxSnow_altitude =  PARA.ensemble.altitude + PARA.snow.relative_maxSnow;
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
        
        water_fluxes_final=zeros(numlabs,numlabs); % Need to be defined before the "if PARA.modules.lateral" in case the variable is not created and used as an inupt in sum_up_output_store
        
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
                        dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s]
                        heat_fluxes = zeros(1, numlabs);
                        PACKAGE_heatExchange.T = T;
                        PACKAGE_heatExchange.cT_grid = GRID.general.cT_grid;
                        PACKAGE_heatExchange.k_cTgrid = k_cTgrid;
                        for j=1:number_of_realizations
                            if j~=index
                                labSend( PACKAGE_heatExchange, j, 1);
                            end
                        end
                        for j=1:number_of_realizations
                            if j~=index
                                PACKAGE_heatExchange_j = labReceive(j, 1);
                                [dE_dt_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index
                                heat_fluxes(j) = nansum( dE_dt_lateral_j .* GRID.general.K_delta );  % in [ W / m^2 ], only for tracking the overall heat fluxes
                                dE_dt_lateral = dE_dt_lateral + dE_dt_lateral_j;    % sum up contributions from all realizations
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
                        water_fluxes= nan(numlabs,numlabs); % in [m/s]
                        % water_fluxes = nan( 1, numlabs ); 
                        PACKAGE_waterExchange.water_table_altitude = PARA.ensemble.water_table_altitude(index);
                        PACKAGE_waterExchange.active_layer_depth_altitude = PARA.ensemble.active_layer_depth_altitude(index);
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
                        
                        water_fluxes_final=nansum(water_fluxes_gather,3);
                        waterflux=nansum(water_fluxes_final(labindex,:));
                        BALANCE.water.dr_lateral = BALANCE.water.dr_lateral + waterflux*1000;
                        waterflux=waterflux+boundary_water_flux;
                        
                        % apply lateral water flux directly (as bulk subsurface flux)
                        [wc, excess_water, lacking_water] = bucketScheme(T, wc, zeros( size(wc) ), GRID, PARA, waterflux);
                        fprintf('\t\t\tExchanged :\t%3.2e m\n\t\t\tBoundary :\t%3.2e m\n\t\t\tTotal :\t%3.2e m\n\t\t\tlacking : %3.2e m\n', waterflux-boundary_water_flux, boundary_water_flux, waterflux, lacking_water)
                        assert( lacking_water < 1e-9, 'CryoGrid3 - lateral exchange - lacking water>0');    % there should be no lacking water as this was checked for
                        if excess_water>0
                            BALANCE.water.dr_lateralExcess=BALANCE.water.dr_lateralExcess + excess_water;            % Added by Léo to have the lateral fluxes in BALANCE
                            GRID.lake.residualWater = GRID.lake.residualWater + excess_water;         % for now: store the excess water, later: treat it according to BC for surface water fluxes % Léo : bot sure why we store this here and not in a BALANCE variable ?
                        end
                        
                        % track water fluxes in BALANCE struct 
                        if strcmp(PARA.ensemble.boundaryCondition(labindex).type,'DarcyReservoir')==1;
                            BALANCE.water.dr_DarcyReservoir = BALANCE.water.dr_DarcyReservoir + boundary_water_flux*1000;
                        end
                        BALANCE.water.dr_water_fluxes_out=BALANCE.water.dr_water_fluxes_out+water_fluxes_final;
                        % for debugging: print water flux per column
                        % fprintf('Water flux to worker %d = %f mm \n', [index, waterflux*1000] );
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
                        mobile_snow = zeros( 1, numlabs );
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
                            BALANCE.water.dr_lateralSnow = BALANCE.water.dr_lateralSnow + my_snow_change ;
                        end
                        
                        snow_fluxes = zeros( 1, numlabs );
                        snow_fluxes(index) = my_snow_change ./ (PARA.technical.syncTimeStep.*3600.*24); %this is only to have comparable output to other fluxes
                        fprintf('\t\t\tSnow flux to worker %d = %f mm SWE \n', [index, my_snow_change*1000] );
                        
                    end
                end
                
                labBarrier();
                PARA = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA) ;
                
                TEMPORARY.syncTime=round((TEMPORARY.syncTime + PARA.technical.syncTimeStep)./PARA.technical.syncTimeStep).*PARA.technical.syncTimeStep;
                fprintf('\t\t\tsync - done\n\n');
            end
        end
        
        
        %------- next time step -----------------------------------------------
        t=t+timestep;
        %---------- sum up + OUTPUT -------------------------------------------
        [TEMPORARY, OUT, BALANCE] = sum_up_output_store(t, T, wc, lwc_cTgrid(GRID.soil.cT_domain), timestep, TEMPORARY, BALANCE, PARA, GRID, SEB, OUT, run_number, snow_fluxes, heat_fluxes);
        
        
    end
    
    % save final state and output at t=endtime
    iSaveOUT(['./runs/' run_number '/' run_number '_output' datestr(t,'yyyy') '.mat'], OUT)
    iSaveState(['./runs/' run_number '/' run_number '_finalState' datestr(t,'yyyy') '.mat'], T, wc, t, SEB, PARA, GRID)
    
end

if number_of_realizations>1
    delete(gcp('nocreate'))
end

disp('Done.');