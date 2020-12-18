function CryoGrid3_xice_mpi_IS(PARA)
% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
% Development by: S. Westermann and M. Langer 2015
% Development of the parallelized version by: L. Martin and J. Nitzbon 2017-2018
% Modified for linear infrastructure by T.Schneider v.D. 2018-2020
% Extended for 2D (fuel storage tank) June 2020

PARA = loadExperimentSetting( PARA ); % load model parameters and settings for experiment  tsvd IS 

numTiles = numlabs;       
SetLoc = PARA.Exp.SetLoc; % specify location (current choices: Prudhoe, Happy Valley, Drew Point)
Scen = PARA.Exp.Scen; % specify scenario (current choices: RCP85, CTR)
Model='GFDL'; %Model='ERA';
%PARA.forcing.filename=[Model,'_',SetLoc,'_',Scen,'_1970-2100'];
PARA.forcing.filename = PARA.Exp.forcing.filename;
% switches for modules
PARA.modules.infiltration=1;    % true if infiltration into unfrozen ground occurs
% PARA.modules.xice=0;          % now specified in loadExperimentSetting.m 
PARA.modules.lateral=1;         % true if adjacent realizations are run (this does not require actual lateral fluxes)  
if(PARA.modules.parallelMode==0); PARA.modules.lateral=0; end  % switch of lateral fluxes for single mode
if PARA.modules.lateral
    % switches for lateral processes
    PARA.modules.exchange_heat = 1;
    PARA.modules.exchange_water = 1;  
    PARA.modules.exchange_snow = 0;
else
    PARA.modules.exchange_heat=0; PARA.modules.exchange_water=0; PARA.modules.exchange_snow=0; 
end
% tsvd IS   specify output directory
%OutDir = '/data/permarisk/CryoGrid3/Runs_paper_final/'
%OutDir = '/data/permarisk/CryoGrid3/Runs_Norilsk/'
OutDir=PARA.Exp.OutDir ;
ExpSet=PARA.Exp.ExpSet;
if(PARA.modules.xice); ExpSet=[ExpSet,'_xice'];end

disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp([' Running Experiment ',ExpSet]) 
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
run_number = [SetLoc,'_',Scen,'_',ExpSet];
%    saveDir = [OutDir,ExpSet,'/',TileSet ]
saveDir = [OutDir,ExpSet,'/',num2str(numTiles),'tiles' ]
if ~exist([saveDir '/Restart'],'dir'); mkdir([saveDir '/Restart']); end %tsvd IS  save workspace for restart
%%%restartFile = [saveDir,'/Restart/restart_',run_number,'_T',num2str(labindex)];
restartDir=[saveDir '/Restart/']; restartFile=['restart_',SetLoc,'_',Scen,'_',ExpSet,'_T',num2str(labindex)];
 
%diary on; diary([ saveDir '/' run_number '_T' num2str(labindex) '_log.txt'])
if(PARA.modules.restart) % run from restart
    load([restartDir,restartFile])
    PARA.technical.endtime=datenum( 2076, 1, 1); % set new endtime of the simulation    
    % PARA.technical.endtime=PARA.technical.new_endtime;  % new_endtime not available anymore after load struct PARA...
    ttt=datevec(t);
    PARA.technical.restartyear=[PARA.technical.restartyear ttt(1)]; % tracking of restart years
else     
    %------ model parameters --------------------------------------------------
%tsvd IS    
    % parameters related to soil
    PARA.soil.albedo=0.2;      %tsvd IS  albedo snow-free surface,  re-defined in get_parallel_variables.m for multiple tiles
    PARA.soil.epsilon=0.97;     % emissvity snow-free surface
    PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
    PARA.soil.rs=50;            % surface resistance against evapotransiration [m^-1] snow-free surface
    PARA.soil.rs_frozen=100;    %tsvd  surface resistance against sublimation [m^-1] snow-free surface
    PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
    PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]
    
    % parameters related to hydrology scheme
%tsvd IS
    PARA.soil.fieldCapacity_MS=0.25;        % water holding capacity of sandy soil (must be consistent with layer porosities...) 15-25% FC https://nrcca.cals.cornell.edu/soil/CA2/CA0212.1-3.php
    PARA.soil.fieldCapacity_Peat=0.5;       % water holding capacity of peat soil   
%    PARA.soil.fieldCapacity_Peat=0.2;       % water holding capacity of peat soil    tsvd IS  see Fig.5.4 FC is between 0.1 and 0.2 for decomposed peat (in Physical Properties of Organic Soils Elon S. Verry, Don H. Boelter, Juhani Päivänen, Dale S. Nichols, Tom Malterer, and Avi Gafni   https://www.nrs.fs.fed.us/pubs/jrnl/2011/nrs_2011_verry_003.pdf
    PARA.soil.fieldCapacity_Gravel=0.1;     % water holding capacity of gravel embankment - Neill & Burn, 2017, table2 https://www.nrcresearchpress.com/doi/pdf/10.1139/as-2016-0036
    PARA.soil.fieldCapacity_Gravel_surface=0.1;  % water holding capacity of gravel surface layer in embankment  
    PARA.soil.evaporationDepth=0.10;        % depth to which evaporation occurs - place on grid cell boundaries   
    PARA.soil.rootDepth=0.20;               % depth affected by transpiration - place on grid cell boundaries
    PARA.soil.ratioET=0.5;                  % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.  % gets overwritten in get_parallel_vars.m
% 
    PARA.soil.externalWaterFlux=0.;        % external water flux / drainage in [m/day]  tsvd IS  gets  re-defined tile specific in get_par_vars.m    
%    PARA.soil.externalWaterFlux=0.01;        % external water flux / drainage in [m/day]
    PARA.soil.drainage=0; % if activated, drainage of outermost tundra tile active and defined according to parameter setting in get_par_vars.m
    PARA.soil.convectiveDomain=[];          % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
    PARA.soil.mobileWaterDomain=[0 10.0];   % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile - numbers are arbitrary
%tsvd IS  If ponding scenario, relative_maxWater now is attributed tile-specific in get_parallel_variables.m   
%    PARA.soil.relative_maxWater=0.05; % loaded in LoadExperimentSetting.m   depth at which a water table will form [m] - above excess water is removed, below it pools up
%    PARA.soil.relative_maxWater_ponding=0.2; % ponding depth of tile specific water table (i.e. a water table can form above the water level defined by relative_maxWater)
    PARA.soil.hydraulic_conductivity = 1e-5;% subsurface saturated hydraulic conductivity assumed for lateral water fluxes [m/s]  Williams and Smith book Fig. 7.10 page 198  for T=0°C, conductivity ~1E-6 m/sec
%tsvd    PARA.soil.infiltration_limit_depth=2.0; % maxiumum depth [m] from the surface to which infiltration occurse
    PARA.soil.infiltration_limit_depth=15.0; % maxiumum depth [m] from the surface to which infiltration occurse    - allow for deeper infiltration
%%%%%%%tsvd IS loadSoilTypes now called after get_parallel_variables.m
    PARA = loadSoilTypes( PARA );           % load the soil types ( silt, sand, water body ) new: gravel
    % parameters related to snow
    PARA.snow.max_albedo=0.85;          % albedo of fresh snow
    PARA.snow.min_albedo=0.5;           % albedo of old snow
    PARA.snow.epsilon=0.99;             % surface emissivity snow
    PARA.snow.z0=5e-4;                  % roughness length surface [m]
    PARA.snow.rs=0.0;                   % surface resistance -> should be 0 for snow
  %  PARA.snow.rho_snow=300;             % density in [kg/m3]
    PARA.snow.tau_1=86400.0;            % time constants of snow albedo change (according to ECMWF reanalysis) [sec]
    PARA.snow.tau_a=0.008;              % [per day]
    PARA.snow.tau_f=0.24;               % [per day]
%tsvd IS  now defined in LoadExperimentSetting.m   relative_maxSnow gets updated in get_par_vars.m to define the value tile_specific    
    %PARA.snow.relative_maxSnow=0.4; 	% maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
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
    PARA.technical.starttime=datenum(1979, 6, 1 );     % starttime of the simulation - if empty start from first value of time series
    PARA.technical.endtime=datenum( 2019, 12, 31, 23, 0, 0); % endtime of the simulation - if empty end at last value of time series
    %PARA.technical.endtime=datenum( 1979, 6, 3);      % endtime of the simulation - if empty end at last value of time series
    PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;       % smallest possible time step in [days] - here 0.1 seconds
    PARA.technical.maxTimestep=300 ./ 3600 ./ 24;       % largest possible time step in [days] - here 300 seconds
    PARA.technical.targetDeltaE=1e5;                    % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
    PARA.technical.outputTimestep= 3./24;           % output time step in [days] - here three hours
%tsvd    PARA.technical.syncTimeStep = 6 ./ 24.0;            % output time step in [days] - here three hours
    PARA.technical.syncTimeStep = 3./24 ;            % output time step in [days] - here three hours  - use same time step as for OUT averaging!
    %PARA.technical.syncTimeStep = PARA.technical.maxTimestep;
    if(numTiles>5); PARA.technical.syncTimeStep = 1./24; end % reduce sync time step to 30min for 24 tile setting (i.e. 1m resolution)!
%    PARA.technical.syncTimeStep = 1 ./ 24.0            % output time step in [days] - here three hours  - use same time step as for OUT averaging!
    PARA.technical.saveDate='01.01.';                   % date of year when output file is written - no effect if "saveInterval" is empty 
    PARA.technical.saveInterval=1;                    % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
    PARA.technical.waterCellSize=0.02;                  % default size of a newly added water cell when water ponds below water table [m]
    
    % subsurface grid
    PARA.technical.subsurfaceGrid = [[0:0.02:4], [4.1:0.1:10], [10.2:0.2:20], [21:1:30], [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]
    %PARA.technical.subsurfaceGrid = [[0:0.02:1], [1.1:0.1:10], [10.2:0.2:20], [21:1:30]]';%, [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]  tsvd

    % parameters related to the specific location (gets updated in get_parallel_variables.m according to ensemble parameter settings)
%tsvd PARA.location.area=1.0;                             % area represented by the run [m^2] (here a dummy value of 1.0 is set which is overwritten for individual tiles)
% attention: initial_altitude gets newly defined in get_parallel_variables for all ensemble members!
%     PARA.location.latitude  = 70.2;         % [deg]  Deadhorse                                                                                                                                               -
%     PARA.location.longitude = -148.46;      % [deg]        
%     PARA.location.initial_altitude=0.0;                % altitude in [m a.s.l.]
    PARA.location.latitude  = 52.39;         % [deg]  Potsdam                                                                                                                                               -
    PARA.location.longitude = 13.06;      % [deg]        
    PARA.location.initial_altitude=35.;                % altitude in [m a.s.l.]
    %     if(PARA.modules.infrastructure==1) % Deadhorse - altitude in [m a.s.l.]
%         PARA.location.initial_altitude=20.0;            
%     elseif(PARA.modules.infrastructure==2)
%         PARA.location.initial_altitude=100.0; 
%     end
    % dynamic auxiliary varaibles (stored in the PARA struct for technical reasons)
    PARA.location.surface_altitude=PARA.location.initial_altitude;          % refers to the surface including water body and snow
    PARA.location.altitude = PARA.location.initial_altitude;                % refers to the terrain surface, including water but excluding snow; used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given
    PARA.location.soil_altitude = PARA.location.initial_altitude;           % refers to the soil surface, excluding water body and snow
    PARA.location.infiltration_altitude = nan;                              % defined at runtime
    PARA.location.water_table_altitude = nan;                               % defined at runtime
    PARA.soil.infiltration_limit_altitude=PARA.location.initial_altitude-PARA.soil.infiltration_limit_depth;    % absolute altitude to which infiltration occurs, updated when ground subsides
    PARA.location.bottomBucketSoilcTIndex = nan; % defined at runtime
    
    % common thresholds 
% tsvd IS set to..... ccc
    PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater;
    %%%PARA.location.absolute_maxWater_altitude = 0.;
    if isempty( PARA.snow.relative_maxSnow )
        PARA.location.absolute_maxSnow_altitude = [];
    else
        PARA.location.absolute_maxSnow_altitude = PARA.location.altitude + PARA.snow.relative_maxSnow;
    end
    
    %initial temperature profile -> first column depth [m] -> second column temperature [degree C]
% tsvd IS
%     if(PARA.modules.infrastructure==1) % Deadhorse
%         %PARA.forcing.filename='GFDL_PrudhoeBay_RCP85_1970-2100';
% %         PARA.Tinitial = [   -2     5    ;...
% %                              0     0    ;...
% %                              2    -2    ;...
% %                              5    -5    ;...
% %                             10    -6    ;...
% %                             25    -6    ;...
% %                            100    -5    ;...
% %                           1100    10.2  ];      % the geothermal gradient for Qgeo=0.05W/m^2 and K=2.746W/Km is about 18.2 K/km
% 
%     PARA.Tinitial = [-2     10    ;...
%                              0     5    ;...
%                              0.5   0.    ;...
%                             10    -7.0    ;...
%                             20    -8.6    ;...  % estimated MAGT at 20m at 1970, borehole data for 1980 Deadhorse (MS draft Walker 2020, Fig. 9)
%                            %100    -5    ;...
%                            600     -1    ;...  % (Lachenbruch JGR 82, borehole Prudhoe Bay)
%                            1100    11.5  ];    % Qgeo=0.05W/m^2 and k=2.0 W/Km ( =(0.3*sqrt(0.57)+0.7*sqrt(3.0))^2, results in ~2.5 K/100m temperature increase
%     end
    
    PARA = loadConstants( PARA );   % load natural constants and thermal properties of soil constituents into the PARA struct
    
    %FORCING data mat-file
    %PARA.forcing.filename='PrudhoeBay_ERAinterim_CCSM4_1979_2094';
    PARA.forcing.rain_fraction=1;   % scaling factor applied to the entire snowfall forcing data
    PARA.forcing.snow_fraction=1;   % scaling factor applied to the entire snowfall forcing data
    PARA.forcing.snow_scaling=1.0;  % scaling factor for incoming snowfall of individual tile, used to emulate lateral snow redistribution
    %---------overwrites variables for each realization--------------------
    % this function must define everything that is realization-specific or dependent of all realizations    
    PARA = get_parallel_variables( PARA );   %cccc only call if numTiles>1....!
    
    if(labindex==1) % display setting
        disp('******************************************************************************************')
        disp('')
        disp('')
        disp('')
        if(PARA.modules.restart==1); disp('RESTART run !!!'); end
        disp(' Experiment Setting:')
        disp(['xice: ',num2str(PARA.modules.xice)])
        disp(['Lateral off/on: ',num2str(PARA.modules.lateral)])
        disp(['Lateral heat: ',num2str(PARA.modules.exchange_heat)])
        disp(['Lateral water: ',num2str(PARA.modules.exchange_water)])
        disp(['embankement height above ground: ',num2str(PARA.IS.EBHag)])
        disp(['ExpSet: ', ExpSet])
        pause on; pause(30)
    end

    %--------------------------------------------------------------------------
    %-----------do not modify from here onwards--------------------------------
    %--------------------------------------------------------------------------
    [FORCING, success]=load_forcing_from_file(PARA); % load FORCING mat-file

    if ~success
        warning('A problem with the Forcing occured.');
    end
%tsvd IS   update forcing...
%% Gravel Road
    if strcmp(PARA.Exp.Case,'GravelRoad') % Gravel Road
            switch string(PARA.IS.TileType)  % remove snow from road and re-distribute to embankment shoulder and toe
                case 'road' % assume that snow is always removed from road
                    FORCING.data.snowfall = zeros(length(FORCING.data.snowfall),1);  % assume that snow is always removed imediately (i.e. prevent snowfall) 
                case 'shoulder'  
                     FORCING.data.snowfall =  FORCING.data.snowfall *4.;
                     if(PARA.IS.RoadOrientation) % use increased SWR for south-facing road
                         FORCING.data.Sin = FORCING.data.Sin_south;
                     end  
                case 'toe'  
                     FORCING.data.snowfall =  FORCING.data.snowfall *4.; 
            end
    end
%% Fuel Tank
    if strcmp(PARA.Exp.Case,'FuelTank')
            switch string(PARA.IS.TileType)  % remove snow from road and re-distribute to embankment shoulder and toe
                case {'pile','tank_bottom'}    
                    FORCING.data.Tair = FORCING.data.TairLowPass0p005;
                    FORCING.data.snowfall = zeros(length(FORCING.data.snowfall),1);
                    FORCING.data.rainfall = zeros(length(FORCING.data.snowfall),1);
                    FORCING.data.wind = FORCING.data.wind/10.; % reduce wind by 1 order of magnitude
                    FORCING.data.Sin = zeros(length(FORCING.data.snowfall),1);
                    %Lout = PARA.surf.epsilon.*sigma.*(T(GRID.air.cT_domain_lb+1)+273.15).^4 + (1-PARA.surf.epsilon).*FORCING.i.Lin;
                    FORCING.data.Lin = PARA.constants.sigma .* (FORCING.data.Tair+273.15).^4;
                %   FORCING.data.q = zeros(length(FORCING.data.snowfall),1);  % leave pressure and specific humidity un-modified   zzz modify q? consistency constrain...?
                case {'shoulder','foundation_base'}    
                    FORCING.data.rainfall = zeros(length(FORCING.data.snowfall),1); % prevent infiltration
                case 'toe'  
                    FORCING.data.rainfall = zeros(length(FORCING.data.snowfall),1); % prevent infiltration
                    FORCING.data.snowfall =  FORCING.data.snowfall *1.5;             % increase snowfall by 50% for snow accumulation at shoulder and toe (retention basin)
            end
    end
%%
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
    %GRID.soil.E_lb = find(PARA.soil.evaporationDepth==GRID.soil.soilGrid(:,1))-1;
    %GRID.soil.T_lb = find(PARA.soil.rootDepth==GRID.soil.soilGrid(:,1))-1;
    
	GRID.soil.water2pool=0; % Leo : cannot find a good initialize function to put it in

    %---- preallocate temporary arrays for capacity and conductivity-----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA); % this is basically the same as "getThermalProperties" during integration, but without interpolation to K grid
    
    %---- energy and water balance initialization -----------------------------
    BALANCE = initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);
    
    %__________________________________________________________________________
    %-------- provide arrays for data storage ---------------------------------
    [t, TEMPORARY] = generateTemporary(T, PARA);
    OUT = generateOUT();
    
    disp('initialization successful');
 %tsvd   iSaveSettings(  [ saveDir '/' run_number '/' run_number '_T' num2str(index) '_settings.mat'] , FORCING, PARA, GRID)
     iSaveSettings(  [ saveDir '/' run_number '_T' num2str(labindex) '_settings.mat'] , FORCING, PARA, GRID)

end % end restart condition
             
    %% ________________________________________________________________________
    % Time Integration Routine                                                I
    while t<PARA.technical.endtime             
        %------ interpolate forcing data to time t ----------------------------
        FORCING = interpolateForcingData(t, FORCING);
        
        %------determine the thermal properties of the model domains ----------
        [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = getThermalPropertiesInfiltration(T, wc, c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid, GRID, PARA);
        
        %------- water and energy balance calculations ------------------------
        BALANCE = updateBALANCE(T, wc, c_cTgrid, lwc_cTgrid, BALANCE, GRID, PARA);
                    assert((sum(isnan(T))==0),['error in T vector AAA for time ',num2str(t)]);        
        
        %------ surface energy balance module ---------------------------------
        %set surface conditions (albedo, roughness length, etc.)
        [PARA, GRID] = surfaceCondition(GRID, PARA, T);
        %calculate the surface energy balance
        [SEB, dwc_dt] = surfaceEnergyBalanceInfiltration(T, wc, FORCING, GRID, PARA, SEB);
        
        %------ soil module  --------------------------------------------------
        %calculate heat conduction
        SEB = heatConduction(T, k_Kgrid, GRID, PARA, SEB);
                    assert((sum(isnan(T))==0),['error in T vector BBB for time ',num2str(t)]);        

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
                    assert((sum(isnan(T))==0),['error in T vector CCC for time ',num2str(t)]);        

        % set grid cells in air to air temperature
        T(GRID.air.cT_domain)=FORCING.i.Tair;
                    assert((sum(isnan(T))==0),['error in T vector DDD for time ',num2str(t)]);        

        %------- water body module --------------------------------------------
        T = mixingWaterBody(T, GRID);
        
        %------- snow cover module --------------------------------------------
        [T, GRID, PARA, SEB, BALANCE] = CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_cTgrid, timestep, BALANCE);
         %assert( ~isnan( GRID.general.K_grid(GRID.snow.cT_domain_ub) ),'error in uppermost snow cell position');        
         assert((sum(isnan(GRID.general.K_grid))==0),'error in GRID.general.K_grid');        
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
%tsvd IS ....            
            %%% PARA.location.absolute_maxWater_altitude = PARA.location.altitude + PARA.soil.relative_maxWater; !tsvd must be soil altitude!
            PARA.location.absolute_maxWater_altitude = PARA.location.soil_altitude + PARA.soil.relative_maxWater;
            %PARA.location.absolute_maxWater_altitude = 0.; %cccc

%tsvd IS   Do not use max soil altitude of ensemble but individual altitudes for IS setting! 
            if isempty( PARA.snow.relative_maxSnow ) 
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
        % all functions called in this block should go into /modules/cryoGridLateral
        % calling PARA.ensemble is only allowed here
    if PARA.modules.lateral
	   if t==TEMPORARY.syncTime %communication between workers
            %tttt fprintf('\n\t\t\tCryoGridLateral: sync - start (Worker %1.0f)\n', labindex);                                
            % update auxiliary variables and common thresholds
            labBarrier();          
            assert((sum(isnan(T))==0),['error1 in T vector for time ',num2str(t)]);        

            [PARA] = updateAuxiliaryVariablesAndCommonThresholds(T, wc, GRID, PARA);

            assert((sum(isnan(T))==0),['error2 in T vector for time ',num2str(t)]);        
                
            % HEAT exchange module
            if PARA.modules.exchange_heat
%tsvd               [ T, TEMPORARY ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid );
                assert((sum(isnan(T))==0),['error3 in T vector for time ',num2str(t)]);        
                [ T, TEMPORARY, BALANCE ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid );
                assert((sum(isnan(T))==0),['error4 in T vector for time ',num2str(t)]);                          
            end
                
            % WATER exchange module
            if PARA.modules.exchange_water
                labBarrier();
                [wc, GRID, BALANCE] = CryoGridLateralWater( PARA, GRID, BALANCE, T, wc);    
            end
                
            % SNOW exchange module
            if PARA.modules.exchange_snow
                PARA  = CryoGridLateralSnow( PARA, GRID );
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
        tt=datevec(t); year=tt(1); %tsvd
        t=t+timestep;
        tt=datevec(t); year_next_ts=tt(1);

        %---------- sum up + OUTPUT -------------------------------------------
        [TEMPORARY, OUT, BALANCE] = sum_up_output_store(t, T, wc, lwc_cTgrid(GRID.soil.cT_domain), timestep, TEMPORARY, BALANCE, PARA, GRID, SEB, OUT, saveDir, run_number);
        if(year~=year_next_ts)
            save([restartDir,restartFile])
        end % save restart file at end of year - gets overwritten for each new year

    end % end while t-loop    
%tsvd  save of outputs done in sum_up_outputstore.m        
fprintf('Done.\n');
end % end function
