add_modules_function;

startDate=datenum( 1979, 6, 1);
endDate=datenum( 1982, 7, 1);

rainFrac=1;
snowFrac=1;
waterTable = 0;
waterTables = [ 0.0, 10.0 ];
maxSnow = 0.4;
snowDens = 200;
extFlux = 2e-3;
fieldCapacity = 0.3;
exices = [ 0.9, 0.8, 0.7, 0.6, 0.5, 0.4];
natPor = 0.4;

jobName = 'TESTRUN';
taskName = sprintf(  [ jobName '_' datestr( startDate, 'yyyymm' ) '-' datestr( endDate, 'yyyymm' ) '_stratSam_rf%d_sf%d_maxSnow%0.1f_snowDens=%0.1f_wt%0.1f_extFlux%0.4f_fc%0.2f' ], ...
                  [ rainFrac, snowFrac, maxSnow, snowDens, ...
                    waterTable, extFlux, fieldCapacity ] );
CryoGrid3_function_spinup(taskName, startDate, endDate, rainFrac, snowFrac, waterTable, maxSnow, snowDens, extFlux, fieldCapacity );