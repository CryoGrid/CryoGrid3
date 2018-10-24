function [t, TEMPORARY] = generateTemporary(T, PARA)

t=PARA.technical.starttime;

TEMPORARY.outputTime=t+PARA.technical.outputTimestep;
TEMPORARY.syncTime=t+PARA.technical.syncTimeStep;

if ~isempty(PARA.technical.saveInterval)
    TEMPORARY.saveTime=datenum(str2num(datestr(t,'yyyy'))+1,  str2num(PARA.technical.saveDate(4:5)), str2num(PARA.technical.saveDate(1:2)))-PARA.technical.outputTimestep;
else
    TEMPORARY.saveTime=PARA.technical.endtime+PARA.technical.outputTimestep;
end

TEMPORARY.Qh_sum=0;
TEMPORARY.Qe_sum=0;
TEMPORARY.Qnet_sum=0;
TEMPORARY.Qg_sum=0;

%for EB checks
TEMPORARY.Qsurf_sum = 0;
TEMPORARY.dE_dt_SEB_sum = 0.*T;
TEMPORARY.dE_dt_cond_sum = 0.*T;

TEMPORARY.dE_soil_sens = 0;
TEMPORARY.dE_soil_lat = 0;
TEMPORARY.dE_soil = 0;
TEMPORARY.dE_snow_sens = 0;
TEMPORARY.dE_snow_lat = 0;
TEMPORARY.dE_snow = 0;


TEMPORARY.timestep_sum=0;

TEMPORARY.T_sum=0.*T;

TEMPORARY.t_last=t;
TEMPORARY.counter=0;

TEMPORARY.waterTableElevation=0;
TEMPORARY.bottomBucketSoilDepth=0;
TEMPORARY.waterTableElevation_sum=0;
TEMPORARY.bottomBucketSoilDepth_sum=0;
TEMPORARY.water2pool_sum=0;
TEMPORARY.water2pool=0;


%TEMPORARY.water_fluxes = zeros( 1, numlabs );     % vector containing accumulated lateral water fluxes per output interval in [m] to the current worker from all other workers
TEMPORARY.snow_flux_lateral = 0 ;      % vector containing accumulated lateral snow fluxes per output interval in [m SWE] to the current worker
TEMPORARY.dE_cell_lateral = zeros( length(T), numlabs );    % matrix containing cell-wise, accumulated lateral heat fluxes in [J/m^3] to the current worker
TEMPORARY.dE_tot_lateral = zeros( 1, numlabs ) ;      % vector containing depth-integrated lateral heat fluxes per output interval in [J/m^2] to the current worker
