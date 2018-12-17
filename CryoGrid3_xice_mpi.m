% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
%
% Development by: S. Westermann and M. Langer 2015
%
% Development of the parallelized version by: L. Martin and J. Nitzbon 2017-2018
%
% -------------------------------------------------------------------------

add_modules;  %adds required modules

number_of_realizations=2;

if number_of_realizations>1 && isempty( gcp('nocreate') )
    parpool(number_of_realizations);
end

spmd
    index=labindex;   % number identifying the process; change this to e.g. 1 for single realization (non-parallel) run
    
    %---------------define input parameters------------------------------------
    % here you provide the ground stratigraphy
    % z     w/i     m       o     type porosity
    % default stratigraphy used in publication:
    PARA.soil.layer_properties=[    0.0     0.55    0.05    0.15    1   0.80    ;...
                                    0.5     0.80    0.05    0.15    1   0.80    ;...
                                    3.0     0.50    0.50    0.00    2   0.50    ;...
                                   10.0     0.03    0.97    0.00    1   0.03    ];
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
    PARA.soil.epsilon=0.97;     % emissvity snow-free surface
    PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
    PARA.soil.rs=50;            % surface resistance against evapotransiration [m^-1] snow-free surface
    PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
    PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]
    
    % parameters related to hydrology scheme
    PARA.soil.fieldCapacity=0.50;           % water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
    PARA.soil.evaporationDepth=0.10;        % depth to which evaporation occurs - place on grid cell boundaries
    PARA.soil.rootDepth=0.20;               % depth affected by transpiration - place on grid cell boundaries
    PARA.soil.ratioET=0.5;                  % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
    PARA.soil.externalWaterFlux=0.0;        % external water flux / drainage in [m/day]
    PARA.soil.convectiveDomain=[];          % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
    PARA.soil.mobileWaterDomain=[0 10.0];   % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
    PARA.soil.relative_maxWater=1.0;        % depth at which a water table will form [m] - above excess water is removed, below it pools up
    PARA.soil.hydraulic_conductivity = 1e-5;% subsurface saturated hydraulic conductivity assumed for lateral water fluxes [m/s]
    PARA.soil.infiltration_limit_depth=2.0; % maxiumum depth [m] from the surface to which infiltration occurse
    PARA = loadSoilTypes( PARA );           % load the soil types ( silt, sand, water body )
    
    % parameters related to snow
    PARA.snow.max_albedo=0.85;          % albedo of fresh snow
    PARA.snow.min_albedo=0.5;           % albedo of old snow
    PARA.snow.epsilon=0.99;             % surface emissivity snow
    PARA.snow.z0=5e-4;                  % roughness length surface [m]
    PARA.snow.rs=0.0;                   % surface resistance -> should be 0 for snow
    PARA.snow.rho_snow=200;             % density in [kg/m3]
    PARA.snow.tau_1=86400.0;            % time constants of snow albedo change (according to ECMWF reanalysis) [sec]
    PARA.snow.tau_a=0.008;              % [per day]
    PARA.snow.tau_f=0.24;               % [per day]
    PARA.snow.relative_maxSnow=[1.0]; 	% maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.snow.extinction=25.0;          % light extinction coefficient of snow [1/m]
    
    % parameters related to water body on top of soil domain
    PARA.water.albedo=0.07;             % albedo water (parameterization after Wayne and Burt (1954) in surfaceCondition.m)
    PARA.water.epsilon=0.99;            % surface emissivity water
    PARA.water.rs=0.0;                  % surface resistance -> should be 0 for water
    PARA.water.z0=5e-4;                 % roughness length surface [m] % value for summer / vegetation    
    PARA.ice.albedo=0.20;               % albedo ice / Lei et al. (2011) shows a range of 0.1 to 0.35
    PARA.ice.epsilon=0.98;              % surface emissivity snow
    PARA.ice.rs=0.0;                    % surface resistance -> should be 0 for ice
    PARA.ice.z0=5e-4;                   % roughness length surface [m] % value for snow
    PARA.ice.extinction=4.5;            % [m^-1] light extinction coefficient of ice / Lei et al. (2011) shows a range of 1 to 5 m^-1
    
    % technical parameters
    PARA.technical.z=2.0;                               % height of input air temperature above ground in [m] - assumed constant even when snow depth increases
    PARA.technical.SWEperCell=0.005;                    % SWE per grid cell in [m] - determines size of snow grid cells
    PARA.technical.maxSWE=0.4;                          % in [m] SWE
    PARA.technical.arraySizeT=5002;                     % number of values in the look-up tables for conductivity and capacity
    PARA.technical.starttime=datenum( 1979, 10, 1 );    % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=datenum( 1979, 12, 31);      % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;       % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;       % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;                    % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3 ./ 24.0;           % output time step in [days] - here three hours
    PARA.technical.syncTimeStep = 6 ./ 24.0;            % output time step in [days] - here three hours
    PARA.technical.saveDate='01.01.';                   % date of year when output file is written - no effect if "saveInterval" is empty
    PARA.technical.saveInterval=[1];                    % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;                  % default size of a newly added water cell when water ponds below water table [m]
    
    % subsurface grid
    %PARA.technical.subsurfaceGrid = [[0:0.02:4], [4.1:0.1:10], [10.2:0.2:20], [21:1:30], [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
    PARA.technical.subsurfaceGrid = [[0:0.02:1], [1.1:0.1:10], [10.2:0.2:20], [21:1:30]]';%, [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]

    % parameters related to the specific location
    PARA.location.area=1.0;                             % area represented by the run [m^2] (here a dummy value of 1.0 is set which is overwritten for individual tiles)
    PARA.location.initial_altitude=20.0;                % altitude in [m a.s.l.]
    
    % dynamic auxiliary varaibles (stored in the PARA struct for technical reasons)
    PARA.location.surface_altitude=PARA.location.initial_altitude;          % refers to the surface including water body and snow
    PARA.location.altitude = PARA.location.initial_altitude;                % refers to the terrain surface, including water but excluding snow; used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given
    PARA.location.soil_altitude = PARA.location.initial_altitude;           % refers to the soil surface, excluding water body and snow
    PARA.location.infiltration_altitude = nan;                              % defined at runtime
    PARA.location.water_table_altitude = nan;                               % defined at runtime
    PARA.soil.infiltration_limit_altitude=PARA.location.initial_altitude-PARA.soil.infiltration_limit_depth;    % absolute altitude to which infiltration occurs, updated when ground subsides
    PARA.location.bottomBucketSoilcTIndex = nan; % defined at runtime
    
    % common thresholds
    PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
    if isempty( PARA.snow.relative_maxSnow )
        PARA.location.absolute_maxSnow_altitude = [];
    else
        PARA.location.absolute_maxSnow_altitude = [ PARA.location.altitude + PARA.snow.relative_maxSnow ];
    end
    
    %initial temperature profile -> first column depth [m] -> second column temperature [degree C]
    PARA.Tinitial = [   -2     5    ;...
                         0     0    ;...
                         2    -2    ;...
                         5    -7    ;...
                        10    -9    ;...
                        25    -9    ;...
                       100    -8    ;...
                      1100    10.2  ];      % the geothermal gradient for Qgeo=0.05W/m^2 and K=2.746W/Km is about 18.2 K/km
    
    PARA = loadConstants( PARA );   % load natural constants and thermal properties of soil constituents into the PARA struct
    
    %FORCING data mat-file
    PARA.forcing.filename='samoylov_ERA_obs_fitted_1979_2014_spinup_extended2044.mat';  %must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
    PARA.forcing.rain_fraction=1;   % scaling factor applied to the entire snowfall forcing data
    PARA.forcing.snow_fraction=1;   % scaling factor applied to the entire snowfall forcing data
    PARA.forcing.snow_scaling=1.0;  % scaling factor for incoming snowfall of individual tile, used to emulate lateral snow redistribution
    
    % switches for modules
    PARA.modules.infiltration=1;    % true if infiltration into unfrozen ground occurs
    PARA.modules.xice=1;            % true if thaw subsicdence is enabled
    PARA.modules.lateral=1;         % true if adjacent realizations are run (this does not require actual lateral fluxes)
    
    if PARA.modules.lateral
        % switches for lateral processes
        PARA.modules.exchange_heat = 1;
        PARA.modules.exchange_water = 0;
        PARA.modules.exchange_snow = 0;
        
        %---------overwrites variables for each realization--------------------
        % this function must define everything that is realization-specific or dependent of all realizations
        PARA = get_parallel_variables( PARA );
    end
    
    % ------make output directory (name depends on parameters) ----------------
    saveDir = './runs';
    run_number = sprintf( [ 'LAT_HEAT_TEST_' datestr( PARA.technical.starttime, 'yyyymm' ) '-' datestr(PARA.technical.endtime, 'yyyymm' )  ] );
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
    [GRID,PARA] = initializeExcessIce(GRID,PARA);
    
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
	GRID.soil.water2pool=0;
    
    %---- preallocate temporary arrays for capacity and conductivity-----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA); % this is basically the same as "getThermalProperties" during integration, but without interpolation to K grid
    
    %---- energy and water balance initialization -----------------------------
    BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);
    
    %__________________________________________________________________________
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    
    disp('initialization successful');
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
        else 
			SEB=surfaceEnergyBalance(T, wc, FORCING, GRID, PARA, SEB);
		end
        
        %------- excess ice module --------------------------------------------
        if PARA.modules.xice && ~PARA.modules.infiltration
            warning( 'energy and water balances are not correct for this combination of modules');
            [GRID, PARA] = excessGroundIce(T, GRID, PARA);
            wc = wc( end-sum(GRID.soil.cT_domain)+1 : end );	% assure wc has correct length
        elseif PARA.modules.xice && PARA.modules.infiltration
            [GRID, PARA, wc, meltwaterGroundIce] = excessGroundIceInfiltration(T, wc, GRID, PARA);
            GRID = updateGRID_excessiceInfiltration(meltwaterGroundIce, GRID);
        end
        
        %------- update Lstar for next time step ------------------------------
        SEB = L_star(FORCING, PARA, SEB);
        
        %------- update auxiliary state variables
        PARA.location.altitude = getAltitude( PARA, GRID );
        PARA.location.surface_altitude = getSurfaceAltitude( PARA, GRID );
        PARA.location.soil_altitude = getSoilAltitude( PARA, GRID );
        [PARA.location.infiltration_altitude, PARA.location.bottomBucketSoilcTIndex] = getInfiltrationAltitude( PARA, GRID, T);
        [PARA.location.water_table_altitude] = getWaterTableAltitudeFC(T, wc, GRID, PARA);
        PARA.soil.infiltration_limit_altitude = PARA.location.soil_altitude - PARA.soil.infiltration_limit_depth;
        
        %------- update threshold variables if no lateral exchange processes occur, otherwise updated at sync time
        if ~PARA.modules.lateral
            PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
            if isempty( PARA.snow.relative_maxSnow ) % Leo : changed from PARA.snow.maxSnow
                PARA.location.absolute_maxSnow_altitude = [];
            else
                PARA.location.absolute_maxSnow_altitude =  PARA.location.altitude + PARA.snow.relative_maxSnow;
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
                                
                % update auxiliary variables and common thresholds
                labBarrier();
                [PARA] = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA);
                
                % HEAT exchange module
                if PARA.modules.exchange_heat
                    [ T, TEMPORARY ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid );
                end
                
                % WATER exchange module
                if PARA.modules.exchange_water
                    [wc, GRID, BALANCE] = CryoGridLateralWater( PARA, GRID, BALANCE, T, wc);    
                end
                
                % SNOW exchange module
                if PARA.modules.exchange_snow
                    PARA  = CryoGridLateralSnow( PARA, GRID );
                end
                
                % lateral EROSION module
                if PARA.modules.exchange_sediment
                    
                    % WRAPPER: CryoGridLateralErosion();
                    
                    sediment_change_tot = zeros(1,numlabs); % in [m] w.r.t. current realization (index)
                    sediment_change_o = zeros(1,numlabs);
                    sediment_change_m = zeros(1,numlabs);


                    PACKAGE_sedimentExchange.K_delta = GRID.general.K_delta( GRID.soil.cT_domain );
                    PACKAGE_sedimentExchange.cT_mineral = GRID.soil.cT_mineral;
                    PACKAGE_sedimentExchange.cT_organic = GRID.soil.cT_organic;

                    
                    
                    for j=1:numlabs
                        if j~=labindex
                            labSend( PACKAGE_sedimentExchange, j, 2);
                        end
                    end
                    for j=1:numlabs
                        if j~=labindex
                            PACKAGE_waterExchange_j = labReceive(j, 2);
                            
                            % WRAPPER: calculte (pair-wise) lateral erosion fluxes
                            
                            if PARA.ensemble.distanceBetweenPoints(labindex,j)>0 
                            
                            
                                soil_surface_index = PARA.ensemble.soil_altitude(labindex);
                                soil_surface_j = PARA.ensemble.soil_altitude(j);

                                surface_index = PARA.ensemble.altitude(labindex);
                                surface_j = PARA.ensemble.altitude(j);

                                cT_mineral_j = PACKAGE_sedimentExchange.cT_mineral ;
                                cT_organic_j = PACKAGE_sedimentExchange.cT_organic ;
                                K_delta_j = PACKAGE_sedimentExchange.K_delta;
                                area_j = PARA.ensemble.area(j);



                                K_delta_index = GRID.general.K_delta(GRID.soil.cT_domain);
                                area_index = PARA.ensemble.area(labindex);

                                K_land = 3e-10; % in [m^2/sec] approx. 0.01 m^2/yr, reference: [ Kessler et al. 2012, JGR ]
                                K_water = 3e-8; % in [m^2/sec] approx 1.0 m^2/yr, reference: [ Kessler et al. 2012, JGR ]              

                                % mixed interface (includes air only)
                                phi_tot = abs( soil_surface_index - soil_surface_j );
                                phi_land = max( soil_surface_index, soil_surface_j ) - min( surface_index, surface_j );
                                phi_land = max( phi_land, 0 );
                                phi_water = min( surface_index, surface_j ) - min( soil_surface_index, soil_surface_j);
                                phi_water = min( phi_water, phi_tot);
                                assert( phi_land + phi_water == phi_tot, 'lateral erosion - water/air interfaces do not match total interface' ) ;
                                K_eff = phi_tot ./ ( phi_land./K_land + phi_water./K_water );       % based on assuming "series junction" of erosion domains --> reciprocal addition of "conductivities"

                                D = PARA.ensemlbe.distanceBetweenPoints(labindex,j);
                                L = PARA.ensemble.thermal_contact_length(labindex,j);

                                % caclulate erosion fluxes
                                q_sed_diff = K_air .* (soil_surface_j - soil_surface_index) ./ D ;  % sediment flux in [m^2/sec]
                                % calculate total sediment volume transported within timestep
                                V_sed_diff = abs( q_sed_diff .* L .* PARA.technical.syncTimeStep .* 24 .* 3600 );  % sediment volume in [m^3] (always positive)                                             


                                if q_sed_diff>0     % gaining sediment

                                    % calculate composition of organic and mineral
                                    V_sed_o = 0;
                                    V_sed_m = 0;
                                    V_sed_remaining = V_sed_diff;
                                    k=1;
                                    while V_sed_remaining > 0
                                        V_temp = min( V_sed_remaining, K_delta_j(k) .* (cT_organic_j(k)+cT_mineral(k)) .* area_j );
                                        V_sed_o = V_sed_o + cT_organic(k) ./ (cT_organic_j(k)+cT_mineral(k)) .* V_temp;
                                        V_sed_m = V_sed_m + cT_mineral(k) ./ (cT_organic_j(k)+cT_mineral(k)) .* V_temp;
                                        V_sed_remaining = V_sed_remaining - V_temp;
                                        k=k+1;
                                    end

                                    sediment_change_tot(j) = V_sed_diff ./ area_index;
                                    sediment_change_m(j) = V_sed_m ./ area_index;
                                    sediment_change_o(j) = V_sed_o ./ area_index;

                                    assert( V_sed_o + V_sed_m == V_sed_diff, 'lateral erosion - organic and mineral do not match total sed change' );



                                elseif q_sed_diff<0     % losing sediment


                                    % calculate composition of organic and mineral
                                    V_sed_o = 0;
                                    V_sed_m = 0;
                                    V_sed_remaining = V_sed_diff;
                                    k=1;
                                    while V_sed_remaining > 0
                                        V_temp = min( V_sed_remaining, K_delta_index(k) .* (GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)) .* area_index );
                                        V_sed_o = V_sed_o + GRID.soil.cT_organic(k) ./ (GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)) .* V_temp;
                                        V_sed_m = V_sed_m + GRID.soil.cT_mineral(k) ./ (GRID.soil.cT_organic(k)+GRID.soil.cT_mineral(k)) .* V_temp;
                                        V_sed_remaining = V_sed_remaining - V_temp;
                                        k=k+1;
                                    end

                                    sediment_change_tot(j) = -V_sed_diff ./ area_index;
                                    sediment_change_m(j) = -V_sed_m ./ area_index;
                                    sediment_change_o(j) = -V_sed_o ./ area_index;

                                    assert( V_sed_o + V_sed_m == V_sed_diff, 'lateral erosion - organic and mineral do not match total sed change' );

                                else % same altitude
                                    sediment_change_tot(j) = 0;
                                    sediment_change_m(j) = 0;
                                    sediment_change_o(j) = 0;

                                end
                                
                            else % not connected
                                sediment_change_tot(j) = 0;
                                sediment_change_m(j) = 0;
                                sediment_change_o(j) = 0;
                            end
                            
                            
                        end
                    end
                    
                    GRID.soil.residual_organic = GRID.soil.residual_organic + sum( sediment_change_o );
                    GRID.soil.residual_mineral = GRID.soil.residual_organic + sum( sediment_change_o );
                    
                    % store fluxes in output (for debugging)
                    
                    
                    % + scale fluxes when >2 realizations

                             
                    % apply erosion fluxes (not implemented so far                 
                    
                    fprintf('\t\t\tsync - lateral erosion\n');
                    fprintf('\t\t\tsync - organic sediment flux to worker %d = %f m \n', [labindex, sum(sediment_change_o) ] );
                    fprintf('\t\t\tsync - mineral sediment flux to worker %d = %f m \n', [labindex, sum(sediment_change_m) ] );

                          
                    
                end
                
                % update auxiliary variables and common thresholds
                labBarrier();
                [PARA] = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA) ;
                
                % determine next sync time
                TEMPORARY.syncTime=round((TEMPORARY.syncTime + PARA.technical.syncTimeStep)./PARA.technical.syncTimeStep).*PARA.technical.syncTimeStep;
                fprintf('\t\t\tsync - done\n');
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
    iPlotAltitudes( [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_altitudes' datestr(t,'yyyy') '.png'], OUT, PARA );
end

if number_of_realizations>1
    delete(gcp('nocreate'))
end

fprintf('Done.\n');
