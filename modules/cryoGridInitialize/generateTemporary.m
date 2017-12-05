function [t, TEMPORARY] = generateTemporary(T, PARA)


t=PARA.technical.starttime;
TEMPORARY.outputTime=t+PARA.technical.outputTimestep;
if ~isempty(PARA.technical.saveInterval)
    TEMPORARY.saveTime=datenum(str2num(datestr(t,'yyyy'))+1,  str2num(PARA.technical.saveDate(4:5)), str2num(PARA.technical.saveDate(1:2)))-PARA.technical.outputTimestep;
else
    TEMPORARY.saveTime=PARA.technical.endtime+PARA.technical.outputTimestep;
end

TEMPORARY.Qh_sum=0;
TEMPORARY.Qe_sum=0;
TEMPORARY.Qnet_sum=0;
TEMPORARY.Qg_sum=0;

TEMPORARY.timestep_sum=0;

TEMPORARY.T_sum=0.*T;

TEMPORARY.t_last=t;
TEMPORARY.counter=0;