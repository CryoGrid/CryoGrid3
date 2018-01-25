% script to excecute multiple independent runs of CryoGrid3 in parallel
% using the job/task batch framwork

add_modules_function;

startDate=datenum( 1979, 6, 1);
endDate=datenum( 1989, 6, 1);

rainFrac=1;
snowFrac=1;
waterTable = 0;
waterTables = [ 0.0, 10.0 ];
maxSnow = 0.4;
maxSnows = [ 0.2, 0.4 ];
snowDens = 200;
snowDensities = [ 200, 250 ];
extFlux = 0;
externalFluxes = [ -2e-3, 0, 2e-3 ];
fieldCapacity = 0.3;
exices = [ 0.9, 0.8, 0.7, 0.6, 0.5, 0.4];
natPor = 0.4;


combinations = {};

i=1;
for ms=maxSnows
    for sd=snowDensities
        for ef=externalFluxes
            combinations{i} = [ms, sd, ef];
            i=i+1;
        end
    end
end

%%

numTasks = length( combinations );

jobName = 'SPINUP';

parallel.defaultClusterProfile('local');
c = parcluster();

job = createJob( c, 'Name', jobName );
disp( [datestr(now) ': created job ' jobName ] );

tasks = {};
for i=1:numTasks
    maxSnow=combinations{i}(1);
    snowDens=combinations{i}(2);
    extFlux = combinations{i}(3);
    
    taskName = sprintf(  [ jobName '_' datestr( startDate, 'yyyymm' ) '-' datestr( endDate, 'yyyymm' ) '_stratSam_rf%d_sf%d_maxSnow%0.1f_snowDens=%0.1f_wt%0.1f_extFlux%0.4f_fc%0.2f' ], ...
                  [ rainFrac, snowFrac, maxSnow, snowDens, ...
                    waterTable, extFlux, fieldCapacity ] );
    tasks{i} = createTask( job , @CryoGrid3_function_spinup, 0 , { taskName, startDate, endDate, rainFrac, snowFrac, waterTable, maxSnow, snowDens, extFlux, fieldCapacity }, 'CaptureDiary', true, 'Name', taskName );
    disp( [ datestr(now) ': created task ' taskName ] );
end

submit(job);
disp( [ datestr(now) ': submitted job ' jobName ] );


% wait(job);
% disp( [ datestr(now) ': finished job ' jobName ] );
% 
% for i=1:numTasks
%     diary( [ './runs/' tasks{i}.Name '/' tasks{i}.Name '_diary.txt' ] );
%     tasks{i}.Diary
%     diary off
%     disp( [datestr(now) ': saved diary of task ' tasks{i}.Name ] );
% end

%delete(job);
%disp( [datestr(now) ': deleted job ' jobName ' - Done.'] );
