% -------------------------------------------------------------------------
% CryoGRID3
% main script for running the model
%
% Developed by: S. Westermann and M. Langer 2015
%
% -------------------------------------------------------------------------

paraFromFile = exist('configFile');     % check if config file passed
add_modules;        % adds required modules
createLogFile=1;    % set to true if the command window output shall be saved

%---------------define input parameters------------------------------------
% here you provide the ground stratigraphy
% z     w/i     m       o     type porosity
PARA.soil.layer_properties=[ 0.0   0.40    0.10    0.15   1   0.75;...   
                             0.15  0.65    0.30    0.05   2   0.65;...   
                             0.9   0.40    0.55    0.05   1   0.40;...    
                             9.0   0.30    0.70    0.00   1   0.30     ];   
% soil stratigraphy
% column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer
% extends until the end of the model domain
% column 2: volumetric water+ice content
% column 3: volumetric mineral content
% column 4: volumetric organic content
% column 5: code for soil type: 1: sand, 2: silt
% column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs

%------ model parameters --------------------------------------------------
% parameters related to surface energy balance and boundary conditions
PARA.soil.albedo=0.2;       % albedo snow-free surface
PARA.soil.albedoPond=0.07;  % albedo of water, used when the uppermost grod cell is 100% water due to modeled thermokarst development
PARA.soil.epsilon=0.97;     % emissvity snow-free surface
PARA.soil.z0=1e-3;          % roughness length [m] snow-free surface
PARA.soil.rs=50;            % surface resistance against evapotransiration [m^-1] snow-free surface
PARA.soil.Qgeo=0.05;        % geothermal heat flux [W/m2]
PARA.soil.kh_bedrock=3.0;   % thermal conductivity of the mineral soil fraction [W/mK]

% parameters related to hydrology scheme
PARA.soil.fieldCapacity=0.3;    %water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
PARA.soil.evaporationDepth=0.1; %depth to which evaporation occurs - place on grid cell boundaries
PARA.soil.rootDepth=0.2;        %depth affected by transpiration - place on grid cell boundaries
PARA.soil.wiltingPoint=0.2;     %point at which transpiration shuts off 
PARA.soil.residualWC=0.05;      %water always remaining in the soil, not accessible to evaporation
PARA.soil.ratioET=0.5;          % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
PARA.soil.externalWaterFlux=0;%-2e-3;  %external water flux / drainage in [m/day]
PARA.soil.convectiveDomain=[];       % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
PARA.soil.mobileWaterDomain=[0 10.0];      % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
PARA.soil.waterTable=0.0;              % depth at which a water table will form [m] - above excess water is removed, below it pools up  

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
PARA.snow.maxSnow= [] ;         % maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
PARA.snow.extinction=25.0;      % light extinction coefficient of snow

% parameters related to the site location
PARA.location.altitude=20.0;    %used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given

% technical parameters
PARA.technical.z=2.0;                       % height of input air temperature above ground in [m] - assumed constant even when snow depth increases
PARA.technical.SWEperCell=0.005;            % SWE per grid cell in [m] - determines size of snow grid cells
PARA.technical.maxSWE=0.4;                  % in [m] SWE
PARA.technical.arraySizeT=5002;             % number of values in the look-up tables for conductivity and capacity
PARA.technical.starttime=datenum('2000.06.01 00:00:00','yyyy.mm.dd HH:MM:SS');       % starttime of the simulation - if empty start from first value of time series
PARA.technical.endtime=datenum('2000.06.15 00:00:00','yyyy.mm.dd HH:MM:SS');         % endtime of the simulation - if empty end at last value of time series
PARA.technical.minTimestep=0.1 ./ 3600 ./ 24;   % smallest possible time step in [days] - here 0.1 seconds
PARA.technical.maxTimestep=300 ./ 3600 ./ 24;   % largest possible time step in [days] - here 300 seconds
PARA.technical.targetDeltaE=1e5;            % maximum energy change of a grid cell between time steps in [J/m3]  %1e5 corresponds to heating of pure water by 0.025 K
PARA.technical.outputTimestep= 3 ./ 24.0 ;          % output time step in [days] - here three hours
PARA.technical.saveDate='01.08.';           % date of year when output file is written - no effect if "saveInterval" is empty
PARA.technical.saveInterval=[];             % interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year
PARA.technical.waterCellSize=0.02;          % default size of a newly added water cell when water ponds below water table [m]
PARA.technical.subsurfaceGrid = [[0:0.02:2], [2.1:0.1:10], [10.2:0.2:20], [21:1:30], [35:5:50], [60:10:100], [200:100:1000]]'; % the subsurface K-grid in [m]

%initial temperature profile -> first column depth [m] -> second column temperature [degree C]
PARA.Tinitial = [-5  10;...
                 0    0;...
                 5    -5;...
                 20    -10;...
                 100   -10;...
                 2000  10];
            
% load natural constants (given in SI units) to PARA.constants
PARA = loadConstants(PARA);

%FORCING data mat-file
PARA.forcing.filename='samoylov_ERA_obs_fitted_1979_2014_spinup.mat';  %must be in subfolder "forcing" and follow the conventions for CryoGrid3 forcing files
PARA.forcing.rain_fraction=1;
PARA.forcing.snow_fraction=1;

% ------update parameter values if config file provided -------------------
% ------changes output directory to name specified in configfile which is the config filename by default
if paraFromFile
    run(configFile);
end

run_number = sprintf('testrun');

% ------make output directory (name depends on parameters) ----------------
mkdir(run_number)

% ------redirect command line output to logfile ---------------------------
if createLogFile
    diary(['./' run_number '/' run_number '_log.mat']);
end


%--------------------------------------------------------------------------
%-----------do not modify from here onwards--------------------------------
%--------------------------------------------------------------------------
[FORCING, success]=load_forcing_from_file(PARA); % load FORCING mat-file

if ~success
    return
end
clear success

PARA = initializeParameters(PARA, FORCING); %set start time, etc.

%----------------create and initialize the grids --------------------------
GRID = makeGrids(PARA);                   %create all grids
GRID = createStratigraphy(PARA,GRID);     %interpolate input stratigraphy to the soil grid

%----- initializie soil thermal properties --------------------------------
GRID = initializeSoilThermalProperties(GRID, PARA);

%------ initializie snow properties----------------------------------------
GRID = initializeSnow(GRID);

%---- initialize the surface energy balance struct ------------------------
SEB = initializeSEB();

%---- initialize temperature profile --------------------------------------
T = inititializeTemperatureProfile_simple(GRID, PARA, FORCING);     

%---- modification for infiltration ---------------------------------------
wc=GRID.soil.cT_water;
GRID.soil.E_lb = find(PARA.soil.evaporationDepth==GRID.soil.soilGrid(:,1))-1;
GRID.soil.T_lb= find(PARA.soil.rootDepth==GRID.soil.soilGrid(:,1))-1;

%---- preallocate temporary arrays for capacity and conductivity-----------
[c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = initializeConductivityCapacity(T,wc, GRID, PARA);

%__________________________________________________________________________
%-------- provide arrays for data storage ---------------------------------
[t, TEMPORARY] = generateTemporary(T, PARA);
OUT = generateOUT();

disp('initialization successful');
save([run_number '/' run_number '_settings.mat'], 'FORCING', 'PARA', 'GRID')

%% ________________________________________________________________________
% Time Integration Routine                                                I
%                                                                         I
%_________________________________________________________________________I
while t<PARA.technical.endtime
        
    %------ interpolate forcing data to time t ----------------------------
    FORCING = interpolateForcingData(t, FORCING);
    
    %------determine the thermal properties of the model domains ----------
    [c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid] = getThermalPropertiesInfiltration(T, wc, c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid, GRID, PARA);

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
    
    %------ determine optimal timestep [days] -----------------------------
    % accounting for min and max timesteps specified, maximum energy change per grid cell and the CFT stability criterion
    timestep = min( [ max( [ min( [ 0.5 * min( GRID.general.K_delta.^2 .* c_cTgrid ./ k_cTgrid ./ (GRID.soil.cT_domain + GRID.snow.cT_domain) ) ./ (24.*3600), ...
                                    PARA.technical.targetDeltaE .* min( abs(GRID.general.K_delta ./ SEB.dE_dt) ) ./ (24.*3600), ...
                                    PARA.technical.maxTimestep ] ), ...
                             PARA.technical.minTimestep ] ), ...
                      TEMPORARY.outputTime-t ] );
    
    %------ update T array ------------------------------------------------
    T = T + SEB.dE_dt./c_cTgrid./GRID.general.K_delta.*timestep.*24.*3600;
    T(GRID.air.cT_domain)=FORCING.i.Tair;
        
    %------- snow cover module --------------------------------------------
    [T, GRID, PARA, SEB] = CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_cTgrid, timestep);
    [GRID, T] = updateGRID_snow(T, GRID, PARA);

    %------- infiltration module-------------------------------------------
    [wc, GRID] = CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING);
      
    %------- update Lstar for next time step ------------------------------
    SEB = L_star(FORCING, PARA, SEB);    
    
    %------- next time step -----------------------------------------------
    t=t+timestep;
    
    %---------- sum up + OUTPUT -------------------------------------------
    sum_up_output_store;    
end
%profile off
save([run_number '/' run_number '_output.mat'], 'OUT')
disp('Done.');