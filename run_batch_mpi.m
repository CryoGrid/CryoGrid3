clear all;
close all;

delete( gcp('nocreate') );

add_modules;

% default setup
SETUP = {};
SETUP.numRealizations = 3;
SETUP.syncTimestep=12./24;

SETUP.startDate = datenum( 2000, 10, 1 );
SETUP.endDate = datenum( 2014, 12, 31 );

SETUP.saveDir = '/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_mpi/runs';


parallel.defaultClusterProfile('local');
c = parcluster();
delete( c.Jobs );

jobName = 'POLYGON-EXPLORE-';
jobs = {};
i=1;

%% default setup
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1e-5;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup with deep ice wedge
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.7;
SETUP.d_xice_T=0.4;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1e-5;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup with shallow ice wedge
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.5;
SETUP.d_xice_T=0.2;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1e-5;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup with different areal fractions
SETUP.f_C=0.4;
SETUP.f_R=0.4;
SETUP.f_T=0.2;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1e-5;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup without external flux
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=0.0;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1e-5;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup with increased evaporation depth
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.10;
SETUP.K=1e-5;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup with lower K [Helbig2013]
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=2e-6;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;

%% setup with higher K [Boike2008]
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1.3e-4;
SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;


%% default setup with long-term run
SETUP.f_C=0.3;
SETUP.f_R=0.6;
SETUP.f_T=0.1;
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T=0.3;
SETUP.extFlux=-0.001;   % in [m/day] this refers to the flux on the landscape-scale and will be adjusted according to the areal fraction
SETUP.d_E=0.02;
SETUP.K=1e-5;

SETUP.startDate = datenum( 1979, 10, 1 );
SETUP.endDate = datenum( 2014, 12, 31 );

SETUP.taskName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) 'fC%0.2f_fR%0.2f_fT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceT%0.2f_extFlux%0.4f_dE%0.2f_K%0.6f' ], ...
    SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T, SETUP.extFlux, SETUP.d_E, SETUP.K ) ;
jobs{i} = batch( @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', 3 );
disp( [datestr(now) ': created task ' SETUP.taskName ] );
i=i+1;


% wait for all jobs to finish
% parfor j=1:i
%     wait( jobs{j} );
%     disp( [ datestr(now) ': finished task number ' num2str(j) ' of ' num2str(i) ] );
%     k=k+1;
% end
