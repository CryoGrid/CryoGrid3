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