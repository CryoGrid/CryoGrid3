function [SETUP] = startFromRunSETUP(startFromRun, suffix)
% Function that generate a setup file to restart a run from a FINAL state
% matfile

if isempty(startFromRun)==0;
    SETUP.flag=1;
    if strcmp(startFromRun(end-3:end),'.mat')==0;
        startFromRun=[startFromRun '.mat'];
    end
    if nargin==1;
        suffix=[];
    end
    
    if ispc
        slash='\';
        k = strfind(startFromRun,'/');
        startFromRun(k)=slash;
    else
        slash='/';
        k = strfind(startFromRun,'\');
        startFromRun(k)=slash;
    end
    
    % load settings to have the number of workers and stuff
    k = strfind(startFromRun,'_realization');
    run_name=startFromRun(1:k-1);
    k= strfind(startFromRun,'_final');
    name_real=startFromRun(1:k-1);
    load(['runs' slash run_name slash name_real '_settings.mat'])
    SETUP.nbreal=length(PARA.ensemble.weight);
    SETUP.run_name_old=run_name;
    SETUP.run_name_new=[run_name '_restart' startFromRun(end-7:end-4) suffix];
    SETUP.restartyear=str2double(startFromRun(end-7:end-4));
    assert(~isnan(SETUP.restartyear),'Error when finding the restart year')
    SETUP.forcingname=PARA.forcing.filename;
    SETUP.slash=slash;
    SETUP.saveDir = ['.' slash 'runs'];
    
    mkdir([ SETUP.saveDir slash SETUP.run_name_new]);
    
else
    
    SETUP.flag=0;
    SETUP.mess='no FINAL state was loaded';
    
end

end

