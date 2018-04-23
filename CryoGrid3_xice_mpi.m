% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
%
% Developed by: S. Westermann and M. Langer 2015
%
% -------------------------------------------------------------------------
    clear all
    close all    
 
%ttt    fileID1 = fopen('log_updateGridInfil_1.txt','w')
%    fileID2 = fopen('log_updateGridInfil_2.txt','w')
    
    par_mode = 1;  % parallel mode off/on
    
if(par_mode==1) 
    delete(gcp('nocreate')) % useful to restart from a crash
end
add_modules;  %adds required modules

%dbstop if error;

number_of_realizations=2;
debug_mode=0   % if set to 1, timestep = timestepMin for debugging (avoid of NaN for timestep calculation)

saveDir = './runs';

if number_of_realizations>1
    parpool(number_of_realizations);
end

%nnn
spmd %zzz use function calls to calls below to enable debugging in par mode!
    index=labindex;   %number identifying the process; change this to e.g. 1 for single realization (non-parallel) run
%nnn    index=1
    
%---------------define input parameters------------------------------------
    % here you provide the ground stratigraphy
    % z     w/i     m       o     type porosity
    
    % default stratigraphy used in publication:
    PARA.soil.layer_properties=[    0.0   0.60    0.10    0.15    1   0.75;...
                                    0.15  0.65    0.3     0.05    2   0.65;...
                                    0.9   0.65    0.3     0.05    1   0.65;...
                                    9.0   0.30    0.70    0.00    1   0.30     ];
                                
%     disp('new stratigraphy........')
%    PARA.soil.layer_properties = [0.0    0.2    0.7    0.00   1   0.30 ;...           
%                                  1.0    0.2    0.7    0.00   1   0.30 ;...
%                                 10.0    0.1    0.8    0.00   1   0.2     ];
    % simple stratigraphy with excess ice used to test water balance:
    % PARA.soil.layer_properties=[ 0.0     0.5    0.5     0.00   1   0.50;...
    %                              0.4     0.8    0.2     0.00   1   0.40;...
    %                              10.0    0.25   0.75    0.00   1   0.25     ];
    % very simply stratigraphy without excess ice used to test energy balance
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
    % parameters related to soil
    PARA.soil.albedo=0.2;       % albedo snow-free surface
    PARA.soil.epsilon=0.97;     % emissvity snow-free surface
    PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
    PARA.soil.rs=50;            % surface resistance against evapotransiration [m^-1] snow-free surface
    PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
    PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]
    
    % parameters related to hydrology scheme
    PARA.soil.fieldCapacity=0.5;    %water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
    PARA.soil.evaporationDepth=0.1; %depth to which evaporation occurs - place on grid cell boundaries
    PARA.soil.rootDepth=0.2;        %depth affected by transpiration - place on grid cell boundaries
    PARA.soil.wiltingPoint=0.2;     %point at which transpiration shuts off
    PARA.soil.residualWC=0.05;      %water always remaining in the soil, not accessible to evaporation
    PARA.soil.ratioET=0.5;          % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
    %tsvd PARA.soil.externalWaterFlux=2e-3;  %external water flux / drainage in [m/day]
    PARA.soil.externalWaterFlux=0.;  %external water flux / drainage in [m/day]
    PARA.soil.convectiveDomain=[];       % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
    PARA.soil.mobileWaterDomain=[0 10.0];      % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
    PARA.soil.relative_maxWater=0.;              % depth at which a water table will form [m] - above excess water is removed, below it pools up
    PARA.soil.hydraulic_conductivity = 1e-5;
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
    PARA.snow.relative_maxSnow= [0.1]; 	%ttt  maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
    PARA.snow.extinction=25.0;      % light extinction coefficient of snow
    %tsvd add lake parameters
    % parameters related to lake
        PARA.water.albedo=0.05;     % albedo water (parameterization after Wayne and Burt (1954) in surfaceCondition.m) 
        PARA.water.epsilon=0.99;    % surface emissivity water
    PARA.water.rs=0.;            % surface resistance -> should be 0 for water
    %tsvd
    PARA.water.z0=5e-4;          % roughness length surface [m] % JAN: value for summer / vegetation
    %tsvd PARA.water.z0=1e-4;         % roughness length surface [m] - gets overridden by value calculated by function flake_roughnessLength.m version Flake
    PARA.water.extinction=1.2;   % light extinction coefficient of water
    PARA.water.depth=0.;
    PARA.water.fetch=20;

    PARA.ice.albedo =0.20;      % albedo ice / Lei et al. (2011) shows a range of 0.1 to 0.35 
    PARA.ice.epsilon=0.98;      % surface emissivity snow
    PARA.ice.rs=0.0;              % surface resistance -> should be 0 for ice
    PARA.ice.z0=5e-4;              % roughness length surface [m] % JAN: value for snow
    PARA.ice.extinction=4.5;    % [m^-1] light extinction coefficient of ice / Lei et al. (2011) shows a range of 1 to 5 m^-1
    
    PARA.technical.z=2.0;                       % height of input air temperature above ground in [m] - assumed constant even when snow depth increases
    PARA.technical.SWEperCell=0.005;            % SWE per grid cell in [m] - determines size of snow grid cells
    PARA.technical.maxSWE=0.4;                  % in [m] SWE
    PARA.technical.arraySizeT=5002;             % number of values in the look-up tables for conductivity and capacity
    PARA.technical.starttime=datenum(1979, 6, 1);       % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=datenum(1980, 12, 31);         % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;   % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;   % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;            % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3 ./ 24.0 ;          % output time step in [days] - here three hours
    PARA.technical.syncTimeStep = 12 ./ 24.0 ;          % output time step in [days] - here three hours
    PARA.technical.saveDate='01.01.';           % date of year when output file is written - no effect if "saveInterval" is empty
    PARA.technical.saveInterval=[];             % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;          % default size of a newly added water cell when water ponds below water table [m]
    
    %default grid used for publications and testing of water balance:
    PARA.technical.subsurfaceGrid = [[0:0.02:4], [4.1:0.1:10], [10.2:0.2:20], [21:1:30], [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
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
	PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
	if isempty( PARA.snow.relative_maxSnow )
    	PARA.location.absolute_maxSnow_altitude = [];
	else
    	PARA.location.absolute_maxSnow_altitude = [ PARA.location.altitude + PARA.snow.relative_maxSnow ];
	end
%tsvd
    PARA.location.latitude  = 72.376;          % [deg]                                                   <                                                                                                        -
    PARA.location.longitude = 126.491;         % [deg]                                                   <                                                                                                        -
    %PARA.location.GridCellSize = 1e6;          % [m^2]                                                   <                                                                                                        -
    %PARA.location.TileFraction = 1e-2;   
    % parameters related to the site location    
    PARA.location.water_table_altitude = nan; % defined at runtime
    
    %initial temperature profile -> first column depth [m] -> second column temperature [degree C]
%     PARA.Tinitial = [ -2     5   ;...
%         0     0   ;...
%         2    -5   ;...
%         10    -10  ;...
%         25    -9   ;...
%         100    -9   ;...
%         2000    10   ];
    PARA.Tinitial = [  -2     5   ;...
                        0     0   ;...
                        2    -2   ;...
                        5    -7   ;...
                        10    -9  ;...
                        25    -9   ;...
                        100    -8   ;...
                        1100    10.2   ];      % the geothermal gradient for Qgeo=0.05W/m² and K=2.746W/Km is about 18.2 K/km 

% load natural constants (given in SI units) to PARA.constants
    PARA = loadConstants( PARA );
    
    %FORCING data mat-file
    PARA.forcing.filename='samoylov_ERA_obs_fitted_1979_2014_spinup.mat';  %must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
    %PARA.forcing.filename='CG3_CCLM_forcing_90_101';
    PARA.forcing.rain_fraction=1;
    PARA.forcing.snow_fraction=1;
    
    % switches for modules
    PARA.modules.infiltration=0;   % true if infiltration into unfrozen ground occurs
    PARA.modules.xice=0;           % true if thaw subsicdence is enabled
	PARA.modules.lateral=1;		   % true if adjacent realizations are run (this does not require actual lateral fluxes)
%tsvd  extended for lateral switched off
    if ~PARA.modules.lateral
        PARA.modules.exchange_heat = 0;
        PARA.modules.exchange_water = 0;
        PARA.modules.exchange_snow = 0;
        PARA.snow.maxSnow= [] ;         % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
        % in par mode this is replaced by PARA.ensemble.immobile_snow_height 
    elseif PARA.modules.lateral
    % switches for lateral processes
        PARA.modules.exchange_heat = 1; %ttt
        PARA.modules.exchange_water = 0; %ttt
        PARA.modules.exchange_snow = 0;  %ttt
        
        %---------overwrites variables for each realization--------------------
		% this function must define everything that is realization-specific or dependent of all realizations
        PARA = get_parallel_variables( PARA );
    end

    disp('Running experiment with xxxx -> indicate switches here')
    % ------make output directory (name depends on parameters) ----------------
%     run_number = sprintf( [ 'Lake-MPI_' datestr( PARA.technical.starttime, 'yyyymm' ) '-' datestr(PARA.technical.endtime, 'yyyymm' ) , ...
%          PARA.modules.exchange_heat, PARA.modules.exchange_water, PARA.modules.exchange_snow, PARA.forcing.rain_fraction, PARA.forcing.snow_fraction, index ] )
    run_number= sprintf( 'LAKE-MPI_xH%d_xW%d_xS%d_infil%d_xice%d_rF%d_sF%d_i%d' , ...
        [ PARA.modules.exchange_heat, PARA.modules.exchange_water, PARA.modules.exchange_snow, ...
          PARA.modules.infiltration, PARA.modules.xice, PARA.forcing.rain_fraction, PARA.forcing.snow_fraction , index ] )
    mkdir([ saveDir '/' run_number]);
    
   %tsvd ------redirect command line output to logfile ---------------------------
   % if createLogFile
%        diary(['./' run_number '/' run_number '_log.mat']);
        diary([ saveDir '/' run_number '/' run_number '_log.mat']);
        % end

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
%tsvd replace by new call lake scheme
%    %---- initialize the water body module ------------------------------------
%    GRID = initializeLAKE(GRID);

%---- initialize the lake scheme structs ---------------------------------- 
    [FLAKE GRID] = initializeLAKE(GRID, PARA);
    %---- initialize temperature profile --------------------------------------
    T = inititializeTemperatureProfile_simple(GRID, PARA, FORCING);
    
    %---- modification for infiltration
    wc=GRID.soil.cT_water;  
    %tsvd wc=GRID.general.cT_water;    
    GRID.soil.E_lb = find(PARA.soil.evaporationDepth==GRID.soil.soilGrid(:,1))-1;
    GRID.soil.T_lb= find(PARA.soil.rootDepth==GRID.soil.soilGrid(:,1))-1;
    
    %---- preallocate temporary arrays for capacity and conductivity-----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA); % this is basically the same as "getThermalProperties" during integration, but without interpolation to K grid
    
    %---- energy and water balance initialization -----------------------------
    BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);
    
    %---- temporary arrays for storage of lateral fluxes --> these could go into a LATERAL struct or TEMPORARY?
%nnn comment out for single mode...
    water_fluxes = zeros( numlabs, 1 );                         % total water flux in [m/s] from each worker to worker index
    snow_fluxes = zeros( numlabs, 1 );                          % total snow flux in [m SWE] per output interval from each worker to worker index
    heat_fluxes = zeros( numlabs, 1);                           % total heat flux in [J] per output interval of all workers to worker index
    dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);    % cell-wise heat flux in [W/m^3]? of all workers to worker index
    
    %__________________________________________________________________________
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    
    disp('initialization successful');
    %%% nnn 
    iSaveSettings(  [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_settings.mat'] , FORCING, PARA, GRID) %tsvd non-interpolated forcing data saved!
    
    %% ________________________________________________________________________
    % Time Integration Routine                                                I
    %                                                                         I
    %_________________________________________________________________________I
    
    while t<PARA.technical.endtime
%tsvd  zzz
    % store old STATE for energy / water balance checks
    T_old = T;
    lwc_old = lwc_cTgrid(GRID.soil.cT_domain);
    wc_old = wc;
    c_cTgrid_old = c_cTgrid;
    Snow_w_old = GRID.snow.Snow_w;

        %------ interpolate forcing data to time t ----------------------------
        FORCING = interpolateForcingData(t, FORCING);
                
        %------determine the thermal properties of the model domains ----------
        [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = getThermalPropertiesInfiltration(T, wc, c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid, GRID, PARA);
        
        %------- water and energy balance calculations ------------------------
        BALANCE = updateBALANCE(T, wc, c_cTgrid, lwc_cTgrid, BALANCE, GRID, PARA);
        
        %------ surface energy balance module ---------------------------------
        %set surface conditions (albedo, roughness length, etc.)
    
%tsvd    [PARA, GRID] = surfaceCondition(GRID, PARA, T); % old implementation
    [PARA, GRID] = surfaceCondition(GRID, PARA, T, t, FORCING, SEB);
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

 %tsvd new timestep calculation to avoid problems with lateral mode switched off
        timestep = min( [ max( [ min( [ 0.5 * nanmin(GRID.general.K_delta.^2 .* c_cTgrid ./ k_cTgrid ./ (GRID.soil.cT_domain + GRID.snow.cT_domain ) ) ./ (24.*3600), ...
                    PARA.technical.targetDeltaE .* nanmin(abs(GRID.general.K_delta ./ SEB.dE_dt ) ) ./ (24.*3600), ...
                    PARA.technical.maxTimestep ] ), ...
                    PARA.technical.minTimestep ] ), ...
                    TEMPORARY.outputTime-t ] );                
        if PARA.modules.lateral
          timestep = min(timestep,TEMPORARY.syncTime-t);
        end
        if(debug_mode)
            timestep =  PARA.technical.minTimestep; % use for debugging to avoid NaN... zzz
        end
        % give a warning when timestep required by CFT criterion is below the minimum timestep specified
        if timestep > 0.5 * min( GRID.general.K_delta.^2 .* c_cTgrid ./ k_cTgrid ./ (GRID.soil.cT_domain + GRID.snow.cT_domain) ) ./ (24.*3600)
            warning( 'numerical stability not guaranteed' );
        end
        
        %------ update T array ------------------------------------------------
        % account for vertical heat fluxes from ground heat flux and heat conduction
        T = T + SEB.dE_dt./c_cTgrid./GRID.general.K_delta.*timestep.*24.*3600;
        % account for lateral heat fluxes
        T = T + dE_dt_lateral./c_cTgrid.*timestep.*24.*3600; % no division by K_delta necessary as dE_dt_lateral in [ J / m^3 / s ]
        % set grid cells in air to air temperature
        T(GRID.air.cT_domain)=FORCING.i.Tair;

        if max((T<-100))==1; disp('dwd'); end

        %------- water body module --------------------------------------------
%tsvd        T = mixingWaterBody(T, GRID);  % zzz check whether this is ok to switch off
        %------- snow cover module --------------------------------------------
        [T, GRID, PARA, SEB, BALANCE] = CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_cTgrid, timestep, BALANCE);
        [GRID, T, BALANCE] = updateGRID_snow(T, GRID, PARA, BALANCE);
        
        %------- infiltration module-------------------------------------------
        if PARA.modules.infiltration
            lateral_water_flux = nansum( water_fluxes );  % lateral fluxes to/from other workers in [m/s]
            [wc, GRID, BALANCE] = CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING, BALANCE, lateral_water_flux);
        end
%tsvd
    %------- water body module --------------------------------------------
    if timestep>0 && ( ~isempty(GRID.lake.water.cT_domain_ub) || ~isempty(GRID.lake.ice.cT_domain_ub) )
        %ice cover
        [T GRID FLAKE] = cryoGridIceCover(GRID, SEB, FLAKE, T, c_cTgrid, timestep);               
        %FLAKE
        if isempty(GRID.lake.ice.cT_domain_ub) && sum(GRID.lake.water.cT_domain)>=3 
            %FLAKE is used if no ice cells exist
            [T GRID FLAKE] = cryoGridFLAKE(FLAKE, GRID, SEB, PARA, T, T_old, timestep, k_Kgrid, SEB.dE_dt, SEB.dE_dt_SEB);
        else
            %updated FLAKE variables also when FLAKE is not used  
            FLAKE.t_bot_n_flk = T(GRID.lake.water.cT_domain_lb) + 273.15;
            %set mixed layerdepth to zero
            FLAKE.h_ml_n_flk = 0;
        end
        GRID = updateGRID_flake(GRID); 
    end 
    %--------- buoyancy calculation in block fields------------------------
    %T = convectionInBlocks(T, c_cTgrid, GRID)  %not available
    % @ ALEX: here the calculations of the excess ice module take place
    % (uncomment this line to "activate" the module). Note that the
    % stratigraphy must be modified such that excess ground ice is present
    % in order to make the module "work"

        %------- excess ice module --------------------------------------------
        if PARA.modules.xice && ~PARA.modules.infiltration
            warning( 'energy and water balances are not correct for this combination of modules');  % zzz ...
            [GRID, PARA] = excessGroundIce(T, GRID, PARA);
%tsvd version Flake [GRID, PARA, wc] = excessGroundIce(T, GRID, PARA);
            % assure wc has correct length
            wc = wc( end-sum(GRID.soil.cT_domain)+1 : end );
%tsvd       wc( end-sum(GRID.soil.cT_domain)+1 : end ); % make vector dims consistent
        elseif PARA.modules.xice && PARA.modules.infiltration
            [GRID, PARA, wc, meltwaterGroundIce] = excessGroundIceInfiltration(T, wc, GRID, PARA);
%tsvd       [GRID, PARA, wc] = excessGroundIceInfiltration(T, wc, GRID, PARA);
            GRID = updateGRID_excessiceInfiltration2(meltwaterGroundIce, GRID);
        end
        
        %------- update Lstar for next time step ------------------------------
        SEB = L_star(FORCING, PARA, SEB);
        
		%------- update auxiliary state variables
		PARA.location.altitude = getAltitude( PARA, GRID );
		PARA.location.surface_altitude = getSurfaceAltitude( PARA, GRID );
 	    PARA.location.active_layer_depth_altitude = getActiveLayerDepthAltitude( PARA, GRID, T );
    	PARA.location.water_table_altitude = getWaterTableAltitude(T, wc, GRID, PARA);

		%------- update threshold variables if no lateral exchange processes occur, otherwise updated at sync time
		if ~PARA.modules.lateral
			PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
            if isempty( PARA.snow.maxSnow )
            %tsvd if isempty( PARA.ensemble.immobile_snow_height )
	   		 	PARA.location.absolute_maxSnow_altitude = [];
			else
				PARA.location.absolute_maxSnow_altitude = [ PARA.ensemble.altitude + PARA.snow.relative_maxSnow ];
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
                disp('CryoGridLateral: sync - start');
                labBarrier(); %common start
                PARA = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA) ; %ddd commenting out (here and below) does not help prevent snow grid crash
                
                % heat exchange module
                if PARA.modules.exchange_heat
                    labBarrier();
                    heat_fluxes = zeros( number_of_realizations, 1);
                    % check preconditions
                    precondition_heatExchange = true; %no specific conditions so far
                    if precondition_heatExchange
                        %WRAPPER
                        disp('sync - exchanging heat');
                        % calculate lateral heat fluxes
                        dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s]
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
                                heat_fluxes(j) = nansum( dE_dt_lateral_j .* GRID.general.K_delta );  % in [J / m^2 s ], only for tracking the overall heat fluxes
                                dE_dt_lateral = dE_dt_lateral + dE_dt_lateral_j;    % sum up contributions from all realizations
                            end
                        end
                    end
                end
                
                % water exchange module
                if PARA.modules.exchange_water
                    labBarrier();
                    water_fluxes = zeros( number_of_realizations, 1 ); % in [m/s]
                    % check preconditions
                    precondition_waterExchange = checkPreconditionWaterExchange( T, GRID );
                    if precondition_waterExchange
                        % WRAPPER
                        disp('sync - exchanging water');
                        % calculate lateral water fluxes
                        PACKAGE_waterExchange.water_table_altitude = PARA.ensemble.water_table_altitude(index);
                        PACKAGE_waterExchange.active_layer_depth_altitude = PARA.ensemble.active_layer_depth_altitude(index);
                        PACKAGE_waterExchange.infiltration_condition = T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub); %jjj
                        for j=1:number_of_realizations
                            if j~=index
                                labSend( PACKAGE_waterExchange, j, 2);
                            end
                        end
                        for j=1:number_of_realizations
                            if j~=index
                                PACKAGE_waterExchange_j = labReceive(j, 2);
                                % JAN: for now: assume DarcyFlux and distribute over sync time step (no check for available water, missing water tracked in BALANCE.dm_lacking)
                                water_fluxes(j) = calculateLateralWaterFluxes( T, PACKAGE_waterExchange_j, GRID, PARA, j);
                            end
                        end
                        % for debugging: print water flux per column
                        waterflux = nansum( water_fluxes.*PARA.technical.syncTimeStep.*24.*3600 ); % in m
                        fprintf('Water flux to worker %d = %f mm \n', [index, waterflux*1000] );
                    end
                    
                end
                
                % snow exchange module
                if PARA.modules.exchange_snow
                    labBarrier();
                    % check preconditions
                    precondition_snowExchange = checkPreconditionSnowExchange( GRID, PARA );
                    if precondition_snowExchange
                        disp('sync - exchanging snow');
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
                            BALANCE.water.dr_lateralSnow = BALANCE.water.dr_lateralSnow + my_snow_change ;
                        end
                        
                        snow_fluxes = zeros( numlabs , 1 );
                        snow_fluxes(index) = my_snow_change ./ (PARA.technical.syncTimeStep.*3600.*24); %this is only to have comparable output to other fluxes
                        fprintf('Snow flux to worker %d = %f mm SWE \n', [index, my_snow_change*1000] );
                        
                    end
                end
                
                labBarrier();
                PARA = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA) ;
                
                TEMPORARY.syncTime=round((TEMPORARY.syncTime + PARA.technical.syncTimeStep)./PARA.technical.syncTimeStep).*PARA.technical.syncTimeStep;
                disp('sync - done');
            end
        end
        
        
        %------- next time step -----------------------------------------------
        t=t+timestep;
        %---------- sum up + OUTPUT -------------------------------------------
      %nnn  
      [TEMPORARY, OUT, BALANCE] = sum_up_output_store(t, T, wc, lwc_cTgrid(GRID.soil.cT_domain), timestep, TEMPORARY, BALANCE, PARA, GRID, SEB, OUT, saveDir, run_number, water_fluxes, snow_fluxes, heat_fluxes);
    end
    
    % save final state and output at t=endtime
%nnn    
iSaveOUT( [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_output' datestr(t,'yyyy') '.mat'], OUT)
%nnn    
iSaveState( [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_finalState' datestr(t,'yyyy') '.mat'], T, wc, t, SEB, PARA, GRID)
%nnn    
%iPlotAltitudes( [ saveDir '/' run_number '/' run_number '_realization' num2str(index) '_altitudes_vs_time_' datestr(t,'yyyy')  '.png'], OUT, PARA );

%nnn 
end

if number_of_realizations>1
    delete(gcp('nocreate'))
end

disp('Done.');

