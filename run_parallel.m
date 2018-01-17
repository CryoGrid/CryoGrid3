% script to excecute multiple indipendent runs of CryoGrid3 in parallel

delete( gcp('nocreate') );

add_modules_function;

number_of_combinations = 12;


startDate=datenum( 1979, 6, 1);
endDate=datenum( 1999, 6, 1);

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

% 
% i=1;
% for wt=waterTables
%     for ex=exices
%         combinations{i} = [wt, ex];
%         i=i+1;
%     end
% end

%start parallel pool
parpool(number_of_combinations, 'SpmdEnabled', false );

parfor i=1:number_of_combinations
    fprintf ( 'Started execution of combination %d ...',  i  );
    %CryoGrid3_function_variableExice( startYear, endYear, rainFrac, snowFrac, combinations{i}(1), maxSnow, snowDens, extFlux, fieldCapacity, combinations{i}(2), natPor);
    CryoGrid3_function_spinup( startDate, endDate, rainFrac, snowFrac, waterTable, combinations{i}(1), combinations{i}(2), combinations{i}(3), fieldCapacity )
    fprintf ( '... finished execution of combination %d.',  i  );
end

delete( gcp('nocreate') );

%parpool(number_of_combinations );
% alternative implementation using spmd
%spmd
   %i=labindex()
   %CryoGrid3_function_spinup( startDate, endDate, rainFrac, snowFrac, waterTable, combinations{i}(1), combinations{i}(2), combinations{i}(3), fieldCapacity )
%end

% alternative implementation using parfeval

% parpool( number_of_combinations );
% % To request multiple evaluations, use a loop.
% for idx = 1:number_of_combinations
%     parfeval(CryoGrid3_function_spinup, 0, startDate, endDate, rainFrac, snowFrac, waterTable, combinations{i}(1), combinations{i}(2), combinations{i}(3), fieldCapacity); % Square size determined by idx
% end
% % Collect the results as they become available.
% % magicResults = cell(1,10);
% % for idx = 1:10
% %   % fetchNext blocks until next results are available.
% %   [completedIdx,value] = fetchNext(f);
% %   magicResults{completedIdx} = value;
% %   fprintf('Got result with index: %d.\n', completedIdx);
% % end

