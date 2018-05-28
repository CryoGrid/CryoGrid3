clear all;
close all;

delete( gcp('nocreate') );

add_modules;

% default setup
SETUP = {};
SETUP.numRealizations = 3;
SETUP.syncTimestep=6./24;

SETUP.startDate = datenum( 1979, 10, 1 );
SETUP.endDate = datenum( 2014, 12, 31);

SETUP.saveDir = './runs';

mkdir( SETUP.saveDir );

jobName = 'THAWSLUMP_';


%% default setup
SETUP.weight_plateau = 100;
SETUP.weight_intermediate = 10;
SETUP.weight_cliff = 1;

% settings for external water reservoir connected to tile "cliff"
SETUP.elevation_Reservoir = -20.0;      % elevation of the reservoir relative to the surface altitude of the "plateau" tile
SETUP.K_Reservoir = 5e-5;               % hydraulic coupling strenth to the reservoir 

SETUP.runName = sprintf( [ jobName datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) ] );



%% excecution in command window

CryoGrid3_xice_mpi(SETUP);


%% alternative: background execution as batch job

parallel.defaultClusterProfile('local');
c = parcluster();

job = batch(c, @CryoGrid3_xice_mpi, 0, { SETUP }, 'CaptureDiary', true, 'Pool', SETUP.numRealizations );
disp( [datestr(now) ': created batch job ' SETUP.runName ] );
