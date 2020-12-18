function run_CryoGrid3_xice_mpi_IS(numTiles) % give number of tiles
% run script to start infrastructure simulations

add_modules;  %adds required modules

PARA.modules.restart = 0; % 0: no restart, 1: restart from restart file (restart file and new end date have to be specified in main script
PARA.technical.restartyear=0; % only for tracking of restart years, do not modify

PARA.modules.parallelMode = 0; % parallel mode off/on

if(numTiles>1)
    PARA.modules.parallelMode = 1; % switch on parallel mode
end

%% infrastructure parameters
PARA.modules.infrastructure=1;  % 1: gravel road setting for Deadhorse
%PARA.IS.EBHag = 2.5;  % embankement height above ground  %todotodo  do switch case ... or load in LoadExpSEtting.m...
%PARA.IS.EBHbg = 1.5;  % embankement height below ground  (total embankment thickness is the sum of EBHag and EBHbg)
PARA.IS.EBHag = 4.0;  % embankement height above ground  %todotodo  do switch case ... or load in LoadExpSEtting.m...
PARA.IS.EBHbg = 0.;  % embankement height below ground  (total embankment thickness is the sum of EBHag and EBHbg)
% further specifications (tile widths, tile tpye, etc. in get_parallel_variables.m
%%

if(PARA.modules.parallelMode==1)
    delete(gcp('nocreate'))
end % must be off in batch mode
if numTiles>1 && isempty( gcp('nocreate') ) % multiple workers
    parpool(numTiles);  % must not be invoked here in batch mode
end

if(PARA.modules.parallelMode==1)
    spmd
    index=labindex;  
   %  diary on; diary([ 'T_',num2str(labindex),'_log.txt'])
    CryoGrid3_xice_mpi_IS(PARA); % call of main program
    end
    delete(gcp('nocreate')) 
else % single worker  - make sure that PARA.modules.lateral set to 0!
   %  diary on; diary([ 'T_',num2str(labindex),'_log.txt'])
    CryoGrid3_xice_mpi_IS(PARA);
end

end % end function
