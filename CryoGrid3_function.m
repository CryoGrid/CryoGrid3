function [] = CryoGrid3_function( runName, startDate, endDate, rainFrac, snowFrac, maxWater, maxSnow, snowDens, extFlux, fieldCapacity,evapDepth, ETversion )

% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
%
% Developed by: S. Westermann and M. Langer 2015
%
% -------------------------------------------------------------------------

paraFromFile = exist('configFile');     % check if config file passed

%add_modules;  %adds required modules

%addpath('./nansuite/')

saveDir = '/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_ETscheme';
mkdir( [ saveDir '/runs' ] );

createLogFile=0;

spinupFile = [];%[ './runs/SPINUP-EXICE_197906-201406_stratSamExice_rf1_sf1_maxSnow0.40_snowDens=200_wt0.0_extFlux0.0000_fc0.30_exice0.60_natPor0.40/SPINUP-EXICE_197906-201406_stratSamExice_rf1_sf1_maxSnow0.40_snowDens=200_wt0.0_extFlux0.0000_fc0.30_exice0.60_natPor0.40_finalState2012.mat' ];

if isempty(spinupFile)
    
    
    %---------------define input parameters------------------------------------
    % here you provide the ground stratigraphy
    % z     w/i     m       o     type porosity
    
    %default used in publication:
    PARA.soil.layer_properties=[ 0.0   0.60    0.10    0.15    1   0.75;...
        0.15  0.65    0.3     0.05    2   0.65;...
        0.9   0.65    0.3     0.05    1   0.65;...
        9.0   0.30    0.70    0.00    1   0.30     ];
    
    %simple stratigraphy with excess ice used to test water balance:
    %PARA.soil.layer_properties=[ 0.0     0.4    0.50    0.00   1   0.50;...
    %                             0.4     0.8    0.20    0.00   1   0.50;...
    %                             10.0    0.25   0.75    0.00   1   0.25     ];
    %very simply stratigraphy without excess ice used to test energy balance
    % soilType = 1;
    % PARA.soil.layer_properties=[ 0.0    0.5    0.5    0.00      soilType  0.5 ;...
    %                             1.0    0.5    0.5    0.00       1   0.5 ;...
    %                             10.0    0.25   0.75   0.00      1   0.25     ];
    % soil stratigraphy
    % column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer
    % extends until the end of the model domain
    % column 2: volumetric water+ice content
    % column 3: volumetric mineral content
    % column 4: volumetric organic content
    % column 5: code for soil type: 1: sand, 2: silt
    % column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs
    
    %------ model parameters --------------------------------------------------
    PARA.soil.albedo=0.2;       % albedo snow-free surface
    PARA.soil.albedoPond=0.07;  % albedo of water, used when the uppermost grod cell is 100% water due to modeled thermokarst development
    PARA.soil.epsilon=0.97;     % emissvity snow-free surface
    PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
    PARA.soil.rs=50;            % surface resistance against evapotransiration [m^-1] snow-free surface
    PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
    PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]
    
    % parameters related to hydrology scheme
    PARA.soil.fieldCapacity=fieldCapacity;    %water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
    PARA.soil.evaporationDepth=evapDepth; %depth to which evaporation occurs - place on grid cell boundaries
    PARA.soil.rootDepth=0.2;        %depth affected by transpiration - place on grid cell boundaries
    PARA.soil.wiltingPoint=0.2;     %point at which transpiration shuts off
    PARA.soil.residualWC=0.05;      %water always remaining in the soil, not accessible to evaporation
    PARA.soil.ratioET=0.5;          % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
    PARA.soil.externalWaterFlux=extFlux;  %external water flux / drainage in [m/day]
    PARA.soil.convectiveDomain=[];       % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
    PARA.soil.mobileWaterDomain=[0 10.0];      % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
    PARA.soil.relative_maxWater=[ maxWater ];              % depth at which a water table will form [m] - above excess water is removed, below it pools up
    PARA = loadSoilTypes( PARA );
    
    % parameters related to snow
    PARA.snow.max_albedo=0.85;      % albedo of fresh snow
    PARA.snow.min_albedo=0.5;       % albedo of old snow
    PARA.snow.epsilon=0.99;         % surface emissivity snow
    PARA.snow.z0=5e-4;              % roughness length surface [m]
    PARA.snow.rs=0.0;               % surface resistance -> should be 0 for snow
    PARA.snow.rho_snow=snowDens;       % density in [kg/m3]
    PARA.snow.tau_1=86400.0;        % time constants of snow albedo change (according to ECMWF reanalysis) [sec]
    PARA.snow.tau_a=0.008;          % [per day]
    PARA.snow.tau_f=0.24;           % [per day]
    PARA.snow.relative_maxSnow= [ maxSnow ]; 	% maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.snow.extinction=25.0;      % light extinction coefficient of snow
    
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
    PARA.technical.starttime=startDate;       % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=endDate;         % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;   % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;   % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;            % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3 ./ 24.0 ;          % output time step in [days] - here three hours
    PARA.technical.saveDate='01.01.';           % date of year when output file is written - no effect if "saveInterval" is empty
    PARA.technical.saveInterval=[1];             % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;          % default size of a newly added water cell when water ponds below water table [m]
    
    %default grid used for publications and testing of water balance:
    PARA.technical.subsurfaceGrid = [[0:0.02:2], [2.1:0.1:10], [10.2:0.2:20], [21:1:30], [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
    %PARA.technical.subsurfaceGrid = [[0:0.02:10], [10.1:0.1:20], [20.2:0.2:30], [31:1:40], [45:5:60], [70:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
    
    PARA.location.area=1.0;
    PARA.location.initial_altitude=20.0;
    % JAN: the following quantities are dynamic and should hence be moved to another struct, e.g. "STATE"
    PARA.location.altitude = PARA.location.initial_altitude; 	% used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given; excluding snow domain
    PARA.location.soil_altitude=PARA.location.initial_altitude;
    PARA.location.surface_altitude=PARA.location.initial_altitude;		% this is dynamic and refers to the surface including snow
    PARA.location.active_layer_depth_altitude = nan; % defined at runtime
    PARA.location.water_table_altitude = nan; % defined at runtime
    % thresholds
    PARA.location.absolute_maxWater_altitude = PARA.location.soil_altitude + PARA.soil.relative_maxWater;
    if isempty( PARA.snow.relative_maxSnow )
        PARA.location.absolute_maxSnow_altitude = [];
    else
        PARA.location.absolute_maxSnow_altitude = [ PARA.location.altitude + PARA.snow.relative_maxSnow ];
    end
    
    %initial temperature profile -> first column depth [m] -> second column temperature [degree C]
    %default:
    PARA.Tinitial = [ -2     5   ;...
        0     0   ;...
        2    -5   ;...
        10    -10  ;...
        25    -9   ;...
        100    -8   ;...
        2000    10   ];
    
    PARA = loadConstants( PARA );
    
    
    %FORCING data mat-file
    PARA.forcing.filename='samoylov_ERA_obs_fitted_1979_2014_spinup.mat';  %must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
    PARA.forcing.rain_fraction=rainFrac;
    PARA.forcing.snow_fraction=snowFrac;
    
    % switches for modules
    PARA.modules.infiltration=1;   % true if infiltration into unfrozen ground occurs
    PARA.modules.xice=1;           % true if thaw subsicdence is enabled
    
    % ------update parameter values if config file provided -------------------
    % ------changes output directory to name specified in configfile which is
    % the config filename by default
    if paraFromFile
        run(configFile);
    end
    
    run_number = runName;
    %     run_number = sprintf( [ 'TESTRUN_' datestr( PARA.technical.starttime, 'yyyymm' ) '-' datestr(PARA.technical.endtime, 'yyyymm' ) '_stratSam_rf%d_sf%d_maxSnow%0.1f_snowDens=%0.1f_wt%0.1f_extFlux%0.4f_fc%0.2f' ], ...
    %                           [ PARA.forcing.rain_fraction, PARA.forcing.snow_fraction, PARA.snow.relative_maxSnow, PARA.snow.rho_snow, ...
    %                           PARA.soil.relative_maxWater, PARA.soil.externalWaterFlux, PARA.soil.fieldCapacity ] );
    
    % ------make output directory (name depends on parameters) ----------------
    mkdir([saveDir '/runs/' run_number])
    
    
    % ------redirect command line output to logfile ---------------------------
    if createLogFile
       % diary([saveDir '/runs/' '/' run_number '_diary.txt']);
    end
    
    
    %--------------------------------------------------------------------------
    %-----------do not modify from here onwards--------------------------------
    %--------------------------------------------------------------------------
    [FORCING, success]=load_forcing_from_file(PARA); % load FORCING mat-file
    
    if ~success
        warning('A problem with the Forcing occured.');
    end
    clear success
    
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
    
    %__________________________________________________________________________
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    
    disp('initialization successful');
    iSaveSettings( [ saveDir '/runs/' run_number '/' run_number '_settings.mat'] , FORCING, PARA, GRID)
    
else %take setting from spinup file
    load(spinupFile); % this loads T, wc, PARA, GRID, SEB from final state of spinup run
    wc(wc<PARA.soil.residualWC)=PARA.soil.residualWC;
    GRID.soil.cT_water = wc;
    GRID = initializeSoilThermalProperties(GRID, PARA);
    
    
    % here one could optionally change the forcing settings
    
    [FORCING, success]=load_forcing_from_file(PARA); % load FORCING mat-file
    if ~success
        warning('A problem with the Forcing occured.');
    end
    clear success
    
    PARA.technical.starttime = datenum(2013, 1, 1);  %take the end time from the spinup run as start time
    PARA.technical.endtime = datenum(2014, 6, 1);
    
    run_number = sprintf( [ 'SPINUP-EXICE_' datestr( PARA.technical.starttime, 'yyyymmdd' ) '-' datestr(PARA.technical.endtime, 'yyyymmdd' ) '_stratSamExice_rf%d_sf%d_maxSnow%0.1f_snowDens=%0.1f_wt%0.1f_extFlux%0.4f_fc%0.2f_exice%0.2f_natPor%0.2f' ], ...
        [ PARA.forcing.rain_fraction, PARA.forcing.snow_fraction, PARA.snow.maxSnow, PARA.snow.rho_snow, ...
        PARA.soil.waterTable, PARA.soil.externalWaterFlux, PARA.soil.fieldCapacity , ...
        PARA.soil.layer_properties(3,2), PARA.soil.layer_properties(3,6) ] );
    
    % ------make output directory (name depends on parameters) ----------------
    mkdir(['./runs/' run_number])
    % necessary initializations
    %---- preallocate temporary arrays for capacity and conductivity-----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA); % this is basically the same as "getThermalProperties" during integration, but without interpolation to K grid
    %---- energy and water balance initialization -----------------------------
    BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    disp('initialization from spinup successful');
    iSaveSettings( [ './runs/' run_number '/' run_number '_settings.mat'] , FORCING, PARA, GRID)
    
    
    
    
    
end
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
    [SEB, dwc_dt] = surfaceEnergyBalanceInfiltration(T, wc, FORCING, GRID, PARA, SEB, ETversion);
    
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
    
    
    % give a warning when timestep required by CFT criterion is below the minimum timestep specified
    %     if timestep > 0.5 * min( GRID.general.K_delta.^2 .* c_cTgrid ./ k_cTgrid ./ (GRID.soil.cT_domain + GRID.snow.cT_domain) ) ./ (24.*3600)
    %         warning( 'numerical stability not guaranteed' );
    %     end
    
    %------ update T array ------------------------------------------------
    T = T + SEB.dE_dt./c_cTgrid./GRID.general.K_delta.*timestep.*24.*3600;
    T(GRID.air.cT_domain)=FORCING.i.Tair;  %set grid cells in air to air temperature
    
    %------- water body module --------------------------------------------
    T = mixingWaterBody(T, GRID);
    
    %------- snow cover module --------------------------------------------
    [T, GRID, PARA, SEB, BALANCE] = CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_cTgrid, timestep, BALANCE);
    
    [GRID, T, BALANCE] = updateGRID_snow(T, GRID, PARA, BALANCE);
    
    %------- infiltration module-------------------------------------------
    if PARA.modules.infiltration
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
    %
    %     % some checks ...
    %     assert( sum( abs( 1 - GRID.soil.cT_actPor - GRID.soil.cT_mineral - GRID.soil.cT_organic )>1e-6 ) == 0, 'CryoGrid3 - actPor has wrong entries' );
    %     assert( sum( sum( GRID.soil.capacity < 0 ) ) == 0, 'CryoGrid3 - negative entry in capacity LUT');
    %     assert( sum( sum( GRID.soil.conductivity < 0 ) ) == 0, 'CryoGrid3 - negative entry in conductivity LUT');
    %
    %
    %------- update Lstar for next time step ------------------------------
    SEB = L_star(FORCING, PARA, SEB);
    
    
    %------- update auxiliary state variables
    PARA.location.altitude = getAltitude( PARA, GRID );
    PARA.location.surface_altitude = getSurfaceAltitude( PARA, GRID );
    PARA.location.soil_altitude = getSoilAltitude( PARA, GRID );
    PARA.location.active_layer_depth_altitude = getActiveLayerDepthAltitude( PARA, GRID, T );
    PARA.location.water_table_altitude = getWaterTableAltitude(T, wc, GRID, PARA);
    %------- update threshold variables
    PARA.location.absolute_maxWater_altitude = PARA.location.soil_altitude + PARA.soil.relative_maxWater;
    if isempty( PARA.snow.relative_maxSnow )
        PARA.location.absolute_maxSnow_altitude = [];
    else
        PARA.location.absolute_maxSnow_altitude = [ PARA.location.altitude + PARA.snow.relative_maxSnow ];
    end
    
    %------- water balance calculations -----------------------------------
    % rainfall
    BALANCE.water.dp_rain = BALANCE.water.dp_rain + FORCING.i.rainfall.*timestep;  %sum up rainfall in [mm] per output interval
    % snowfall
    BALANCE.water.dp_snow = BALANCE.water.dp_snow + FORCING.i.snowfall.*timestep; %sum up snowfall in [mm] SWE per output interval
    
    %------- next time step -----------------------------------------------
    t=t+timestep;
    %---------- sum up + OUTPUT -------------------------------------------
    [TEMPORARY, OUT, BALANCE] = sum_up_output_store(t, T, wc, lwc_cTgrid(GRID.soil.cT_domain), timestep, TEMPORARY, BALANCE, PARA, GRID, SEB, OUT, run_number, saveDir);
    
end

% save final state and output at t=endtime
iSaveOUT([saveDir '/runs/' run_number '/' run_number '_output' datestr(t,'yyyy') '.mat'], OUT)
iSaveState([saveDir '/runs/' run_number '/' run_number '_finalState' datestr(t,'yyyy') '.mat'], T, wc, t, SEB, PARA, GRID)


disp('Done.');
if createLogFile
    diary off;
end
