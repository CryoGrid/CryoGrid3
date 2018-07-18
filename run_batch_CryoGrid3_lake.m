% run CryoGrid3 in batch modus
clear all; close all; delete( gcp('nocreate') );
add_modules;                    % needs to be called in main.m ... (on cluster)
SETUP = {};                        % auxiliary struct containing parameters which are often varied. this struct is passed as an argument to the main script.

parallel.defaultClusterProfile('local'); c = parcluster(); delete( c.Jobs ); % jjj ?
jobName = 'Lateral_Heat_Lake';
jobs = {};                                                                   

%% shallow lake (1-3: SSW, 4-6:MSW)
SETUP.LD = 1;    % lake depth in meter
SETUP.Tini   =  [ -5     15  15;...    % SSW setting
                   0     10    6;...   % surface
                   1    -3    6;...    % lake bottom
                   2    -4    0;...    % talik depth
                  10    -7   -5;...
                  20   -10   -10;...
                 100   -10   -10;...  
                2000    10   10];
% SSW
    SETUP.LR = 10;    % lake radius in meter
% 1) LD 1m  LR 10m  LC=0.1		
i=1;  
SETUP.LF= 0.1;        % lake landscape coverage 
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 2)  LD 1m  LR 10m  LC=0.25				
i=2;    
SETUP.LF= 0.25;        
SETUP.endtime = datenum(2099, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 3)  LD 1m  LR 10m  LC=0.5			
i=3;    
SETUP.LF= 0.5;       
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );

% MSW
    SETUP.LR = 100;    % lake radius in meter
% 4)  LD 1m  LR 100m  LC=0.1				
i=4;    % counter for the jobs
SETUP.LF= 0.1;        % lake landscape coverage 
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 5) LD 1m  LR 100m  LC=0.25		
i=5;    
SETUP.LF= 0.25;        
SETUP.endtime = datenum(2099, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 6) LD 1m  LR 100m  LC=0.5	
i=6;    
SETUP.LF= 0.5;       
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );

%% deep lake  (7-9: SSW, 10-12: MSW)
SETUP.LD = 5;    % lake depth in meter
SETUP.Tini = [  -5     15  15;...  % MSW setting for summer
           0     10   5;...  % soil surface
           5    -3    5;...  % lake bottom
           8    -6    0;...  % talik depth
          20   -10   -9;...
         100   -10   -10;...  
        2000    10   10];     

% SSW
    SETUP.LR = 10;    % lake depth in meter
% 7) LD 5m  LR 10m  LC=0.1		
i=7;    
SETUP.LF= 0.1;        % lake landscape coverage 
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 8) LD 5m  LR 10m  LC=0.25		
i=8;    
SETUP.LF= 0.25;        
SETUP.endtime = datenum(2099, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
%9) LD 5m  LR 10m  LC=0.5	
i=9;    
SETUP.LF= 0.5;       
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );

% MSW
    SETUP.LR = 100;    % lake depth in meter
% 10) LD 5m  LR 100m  LC=0.1		
i=10;    % counter for the jobs
SETUP.LF= 0.1;        % lake landscape coverage 
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 11)  LD 5m  LR 100m  LC=0.25
i=11;    
SETUP.LF= 0.25;        
SETUP.endtime = datenum(2099, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
% 12) LD 5m  LR 100m  LC=0.5	
i=12;    
SETUP.LF= 0.5;       
SETUP.endtime = datenum(2015, 12, 31); 
jobs{i} = batch( @CryoGrid3_lake_batch, 0, { SETUP }, 'CaptureDiary',true, 'Pool', 2 );
disp( [datestr(now) ': created job ',num2str(i) ] );
%%

% delete(job)
% clear job

% % to watch the status use c.Jobs(i)(.Tasks) or jobs{i}(.Tasks)
% % j = batch('CryoGrid3_xice_mpi','Pool',3,'CaptureDiary',true);
% %  wait(j);   % Wait for the job to finish
% %  diary(j)   % Display the diary
% %  load(j)    % Load job workspace data into client workspace
