% script to excecute multiple independent runs of CryoGrid3 in parallel
% using the job/task batch framework

add_modules_function;

startDate=datenum( 1985, 10, 1);
endDate=datenum( 1986, 10, 1);

rainFrac=1;
snowFrac=1;
maxSnow = 1.0;
snowDens = 225;
extFlux = 0;
externalFluxes = [  -2e-3, 0, 2e-3 ];
fieldCapacity = 0.5;
fieldCapacities = [0.3, 0.5];
exices = [ 0.9, 0.8, 0.7, 0.6, 0.5, 0.4];
natPor = 0.4;
maxWater = 0.;
maxWaters = [ 0.5 ];

saturations=[ 0., 0.2, 0.4, 0.6, 0.8, 1.0 ];


% combinations = {};
% 
% i=1;
% for ms=maxSnows
%     for sd=snowDensities
%         for ef=externalFluxes
%             combinations{i} = [ms, sd, ef];
%             i=i+1;
%         end
%     end
% end

combinations = {};

i=1;
for sat=saturations
    for ex=externalFluxes
        combinations{i}= [sat, ex];
        i=i+1;
    end
end

%%

numTasks = length( combinations );

jobName = 'TESTRUN';

parallel.defaultClusterProfile('local');
c = parcluster();

job = createJob( c, 'Name', jobName );
disp( [datestr(now) ': created job ' jobName ] );

tasks = {};
for i=1:numTasks
    saturation=combinations{i}(1);
    extFlux=combinations{i}(2);
    
    taskName = sprintf(  [ jobName '_' datestr( startDate, 'yyyymm' ) '-' datestr( endDate, 'yyyymm' ) '_rf%d_sf%d_maxSnow%0.1f_snowDens=%0.1f_maxWater%0.1f_extFlux%0.4f_fc%0.2f_' ], ...
                  [ rainFrac, snowFrac, maxSnow, snowDens, ...
                    maxWater, extFlux, fieldCapacity ] );
    tasks{i} = createTask( job , @CryoGrid3_function, 0 , { taskName, startDate, endDate, rainFrac, snowFrac, maxWater, maxSnow, snowDens, extFlux, fieldCapacity, evapDepth, ETversion, saturation }, 'CaptureDiary', true, 'Name', taskName );
    disp( [ datestr(now) ': created task ' taskName ] );
end

% for i=1:numTasks
%     exice=combinations{i}(1);
%     waterTable=combinations{i}(2);
%     
%     taskName = sprintf(  [ jobName '_' datestr( startDate, 'yyyymm' ) '-' datestr( endDate, 'yyyymm' ) '_stratSamExice_rf%d_sf%d_maxSnow%0.2f_snowDens=%d_wt%0.1f_extFlux%0.4f_fc%0.2f_exice%0.2f_natPor%0.2f' ], ...
%                    [ rainFrac, snowFrac, maxSnow, snowDens, waterTable, extFlux, fieldCapacity, exice, natPor ] );
% 
%     tasks{i} = createTask( job, @CryoGrid3_function_variableExice, 0, ...
%                             { taskName, startDate, endDate, rainFrac, snowFrac, waterTable, maxSnow, snowDens, extFlux, fieldCapacity, exice, natPor }, ...
%                             'CaptureDiary', true, 'Name', taskName);
%     disp( [ datestr(now) ': created task ' taskName ] );    
% end

submit(job);
disp( [ datestr(now) ': submitted job ' jobName ] );


%wait(job);
%disp( [ datestr(now) ': finished job ' jobName ] );
% 
% for i=1:numTasks
%     diary( [ './runs/' tasks{i}.Name '/' tasks{i}.Name '_diary.txt' ] );
%     tasks{i}.Diary
%     diary off
%     disp( [datestr(now) ': saved diary of task ' tasks{i}.Name ] );
% end

%delete(job);
%disp( [datestr(now) ': deleted job ' jobName ' - Done.'] );
