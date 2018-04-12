clear all;
close all;

delete( gcp('nocreate') );

add_modules;

% default setup
SETUP = {};
SETUP.numRealizations = 3;
SETUP.syncTimestep=12./24;

SETUP.startDate = datenum( 1979, 10, 1 );
SETUP.endDate = datenum( 2014, 10, 5 );

SETUP.saveDir = '/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_mpi_polygon/runs';

mkdir( SETUP.saveDir );

parallel.defaultClusterProfile('local');
c = parcluster();
delete( c.Jobs );

jobName = 'POLYGON_';
jobs = {};
i=1;

%% default setup
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.2;
SETUP.e_Reservoir = 0;
SETUP.K_Reservoir = 1e-5;
SETUP.K=1e-5;
SETUP.relMaxWater = 0.;
SETUP.heatExchangeAltitudeFactor = 1.;

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_maxWater%0.2f_xHeatFactor%d' ], ...
    SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.relMaxWater, SETUP.heatExchangeAltitudeFactor) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created task ' SETUP.runName ] );
i=i+1;

% setup with deep ice wedge
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.7;
SETUP.d_xice_T=0.3;
SETUP.e_Reservoir = 0;
SETUP.K_Reservoir = 1e-5;
SETUP.K=1e-5;
SETUP.relMaxWater = 0.;
SETUP.heatExchangeAltitudeFactor = 1.;

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_maxWater%0.2f_xHeatFactor%d' ], ...
    SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.relMaxWater, SETUP.heatExchangeAltitudeFactor) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created task ' SETUP.runName ] );
i=i+1;

% setup with shallow ice wedge
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.55;
SETUP.d_xice_T=0.15;
SETUP.e_Reservoir = 0;
SETUP.K_Reservoir = 1e-5;
SETUP.K=1e-5;
SETUP.relMaxWater = 0.;
SETUP.heatExchangeAltitudeFactor = 1.;

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_maxWater%0.2f_xHeatFactor%d' ], ...
    SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.relMaxWater, SETUP.heatExchangeAltitudeFactor) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created task ' SETUP.runName ] );
i=i+1;

% setup with low reservoir ice wedge
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.2;
SETUP.e_Reservoir = -0.2;
SETUP.K_Reservoir = 1e-5;
SETUP.K=1e-5;
SETUP.relMaxWater = 0.;
SETUP.heatExchangeAltitudeFactor = 1.;

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_maxWater%0.2f_xHeatFactor%d' ], ...
    SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.relMaxWater, SETUP.heatExchangeAltitudeFactor) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created task ' SETUP.runName ] );
i=i+1;

% setup with high reservoir ice wedge
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.2;
SETUP.e_Reservoir = 0.2;
SETUP.K_Reservoir = 1e-5;
SETUP.K=1e-5;
SETUP.relMaxWater = 0.;
SETUP.heatExchangeAltitudeFactor = 1.;

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_maxWater%0.2f_xHeatFactor%d' ], ...
    SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.relMaxWater, SETUP.heatExchangeAltitudeFactor) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created task ' SETUP.runName ] );
i=i+1;

%setup with heat exchange to top
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.2;
SETUP.e_Reservoir = 0;
SETUP.K_Reservoir = 1e-5;
SETUP.K=1e-5;
SETUP.relMaxWater = 0.;
SETUP.heatExchangeAltitudeFactor = 0.;

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_maxWater%0.2f_xHeatFactor%d' ], ...
    SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.relMaxWater, SETUP.heatExchangeAltitudeFactor) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created task ' SETUP.runName ] );


%wait for all jobs to finish
parfor j=1:i
    wait( jobs{j} );
    disp( [ datestr(now) ': finished task number ' num2str(j) ' of ' num2str(i) ] );
end
