% author: jnitzbon
% description: this script is used to specify individual setups for
% parallel runs of CryoGrid3. it is called by the script "submit_SLURM.sh"
% which is used to submit MATLAB jobs to the SLURM queue
clear all;
close all;

%% generate SETUP struct

SETUP = {};

%startYear=2011;

% parameters
SETUP.numRealizations = 3;
SETUP.syncTimestep=6./24;
SETUP.startDate = datenum( 1979, 10, 1 );
SETUP.endDate = datenum( 2099, 12, 31);
SETUP.xH=1;
SETUP.xW=1;
SETUP.xS=1;
SETUP.xice=1;

SETUP.fieldCapacity = 0.50;

SETUP.e_Reservoir = -1.0;%-1.0;
SETUP.snowDens = 200;%200..250
SETUP.e_R = 0.4;%0.2..0.4
SETUP.e_T = SETUP.e_R-0.1;

SETUP.boundaryCondition_T = 'DarcyReservoir';

% areal fractions
SETUP.f_C = 0.3;
SETUP.f_T = 0.1;
SETUP.f_R = 1.0-SETUP.f_T-SETUP.f_C;

% excess ice layers
SETUP.d_xice_C=0.9;
SETUP.d_xice_R=0.6;
SETUP.d_xice_T1=0.3;
SETUP.d_xice_T2=SETUP.d_xice_R-abs(SETUP.e_R-SETUP.e_T); % intermediate layer with reduced excess ice between T1 and T2, # T2 should be at the same altitude as d_xice_R

% d_xice_DELTA=0.0;
% SETUP.d_xice_C=SETUP.d_xice_C+d_xice_DELTA;
% SETUP.d_xice_R=SETUP.d_xice_R+d_xice_DELTA;
% SETUP.d_xice_T1=SETUP.d_xice_T1+d_xice_DELTA;
% SETUP.d_xice_T2=SETUP.d_xice_T2+d_xice_DELTA;

SETUP.K_Reservoir = 5e-5;
SETUP.K=1e-5;
SETUP.natPor = 0.55;

SETUP.scenario='rcp45';

SETUP.forcingFile = ['Samoylov_' SETUP.scenario '_1901_2300_CryoGrid_windModified.mat'];

% output directory
SETUP.saveDir = '/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_mpi_polygon/runs';

% compose runname
%SETUP.runName = sprintf( [ 'VALIDATION_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_xice%d_xH%d_xW%d_xS%d_fC%0.1f_fR%0.1f_fT%0.1f_eR%0.2f_eT%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_fc%0.2f_snowDens%d' ], ...
 %   SETUP.xice, SETUP.xH, SETUP.xW, SETUP.xS, ...
  %  SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.e_R, SETUP.e_T, ...
   % SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.fieldCapacity, SETUP.snowDens) ;
%SETUP.runName = sprintf( [ 'LONGTERM_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_xice%d_xH%d_xW%d_xS%d_fC%0.1f_fR%0.1f_fT%0.1f_eR%0.2f_eT%0.2f_dxiceC%0.2f_dxiceR%0.2f_dxiceTa%0.2f_dxiceTb%0.2f_K%0.1e_KRes%0.1e_eRes%0.2f_snowDens%d' ], ...
 %    SETUP.xice, SETUP.xH, SETUP.xW, SETUP.xS, ...
  %   SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.e_R, SETUP.e_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T1, SETUP.d_xice_T2, ...
   %  SETUP.K, SETUP.K_Reservoir, SETUP.e_Reservoir, SETUP.snowDens) ;
SETUP.runName = sprintf( [ 'SCENARIO_' SETUP.scenario '_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_xice%d_xH%d_xW%d_xS%d_%s_eRes%0.2f_fC%0.1f_fR%0.1f_fT%0.1f_dxiceC%0.2f_dxiceR%0.2f_dxiceTa%0.2f_dxiceTb%0.2f' ], ...
     SETUP.xice, SETUP.xH, SETUP.xW, SETUP.xS, SETUP.boundaryCondition_T, SETUP.e_Reservoir,...
     SETUP.f_C, SETUP.f_R, SETUP.f_T, SETUP.d_xice_C, SETUP.d_xice_R, SETUP.d_xice_T1, SETUP.d_xice_T2 ) ;
[~, SETUP.git_commit_hash] = system('git rev-parse HEAD');

% create ouput directory
mkdir( [ SETUP.saveDir ] );
mkdir( [ SETUP.saveDir '/' SETUP.runName ] );

% save setup file
save( [ SETUP.saveDir '/' SETUP.runName '/' SETUP.runName '_setup.mat' ], 'SETUP'  );

% create temporary files for output directory and run name
fid = fopen('SAVEDIR.temp','wt');
fprintf(fid, [ SETUP.saveDir ] );
fclose(fid);

fid = fopen('RUN.temp','wt');
fprintf(fid, [ SETUP.runName ] );
fclose(fid);
