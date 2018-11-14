% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
%
% Development by: S. Westermann and M. Langer 2015
%
% Development of the parallelized version by: L. Martin and J. Nitzbon 2017-2018
%
% -------------------------------------------------------------------------

clear
clc
close all
delete(gcp('nocreate')) % useful to restart from a crash

add_modules;  %adds required modules

number_of_realizations=5;

if number_of_realizations>1 && isempty( gcp('nocreate') )
    parpool(number_of_realizations);
end

% Name, Forcing and diary
run_number='181114_50y';
forcingname='Suossjavri_WRF_Norstore_adapted50yr.mat';
diary(['./runs/' run_number '_log.txt'])

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
    PARA.soil.fieldCapacity=0.55;           % water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
    PARA.soil.evaporationDepth=0.05;        % depth to which evaporation occurs - place on grid cell boundaries
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
    PARA.snow.relative_maxSnow=1.0; 	% maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
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
    PARA.technical.starttime=[];                        % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=[];                          % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;       % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;       % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;                    % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3 ./ 24.0;           % output time step in [days] - here three hours
    PARA.technical.syncTimeStep = 6 ./ 24.0;            % output time step in [days] - here three hours
    PARA.technical.saveDate='01.08.';                   % date of year when output file is written - no effect if "saveInterval" is empty
    PARA.technical.saveInterval=1;                      % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;                  % default size of a newly added water cell when water ponds below water table [m]
    
    % subsurface grid
    PARA.technical.subsurfaceGrid = [0:0.05:2 2.1:0.1:5 5.2:0.2:20 21:1:30 35:5:50 60:10:100 200:100:1000]'; % the subsurface K-grid in [m]
    
    % parameters related to the specific location
    PARA.location.area=1.0;                             % area represented by the run [m^2] (here a dummy value of 1.0 is set which is overwritten for individual tiles)
    PARA.location.initial_altitude=300;                % altitude in [m a.s.l.]
    
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
    PARA.Tinitial = [-5     10   ;...
                      0      5   ;...
                      0.1    2   ;...
                      0.5    0   ;...
                      1     -1   ;...
                      2     -2   ;...
                     10      1   ;...
                     30      2   ;...
                    500      4   ;...
                   1000     10   ];     % the geothermal gradient for Qgeo=0.05W/m^2 and K=2.746W/Km is about 18.2 K/km
    
    PARA = loadConstants( PARA );   % load natural constants and thermal properties of soil constituents into the PARA struct
    
    %FORCING data mat-file
    PARA.forcing.filename=forcingname;  %must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
    PARA.forcing.rain_fraction=1;
    PARA.forcing.snow_fraction=0.2;
    
    % switches for modules
    PARA.modules.infiltration=1;    % true if infiltration into unfrozen ground occurs
    PARA.modules.xice=1;            % true if thaw subsicdence is enabled
    PARA.modules.lateral=1;         % true if adjacent realizations are run (this does not require actual lateral fluxes)
    
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
    saveDir = './runs';
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
    GRID.soil.E_lb = find(PARA.soil.evaporationDepth==GRID.soil.soilGrid(:,1))-1;
    GRID.soil.T_lb = find(PARA.soil.rootDepth==GRID.soil.soilGrid(:,1))-1;
    GRID.soil.water2pool=0; % Leo : cannot find a good initialize function to put it in
    GRID.soil.flag.dry2permafrost=0; % Flag to keep track of the existence of a water table above the permafrost 
    GRID.soil.flag.noWTnoPF=0;
    
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
        [PARA.location.water_table_altitude, GRID.soil.flag] = getWaterTableAltitudeFC(T, wc, GRID, PARA);
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
                    [T, GRID, BALANCE, TEMPORARY] = CryoGridLateralSnow( PARA, GRID, BALANCE, TEMPORARY, FORCING, T);
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

fprintf('Done with %s\n', run_number);
diary off