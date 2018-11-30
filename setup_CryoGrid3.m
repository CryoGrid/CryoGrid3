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
SETUP.startDate = datenum( 1949, 10, 1 );
SETUP.endDate = datenum( 2099, 12, 31);
SETUP.xH=1;
SETUP.xW=1;
SETUP.xS=1;
SETUP.xice=0;

SETUP.fieldCapacity = 0.50; % 0.40
SETUP.relMaxSnow = 0.4; % 1.0
SETUP.snowDens = 200;%200..250
SETUP.boundaryCondition_T = 'DarcyReservoir';
SETUP.e_Reservoir = 0.0;%-1.0;

% areal fractions
SETUP.f_C = 0.3; % 0.5
SETUP.f_T = 0.1;
SETUP.f_R = 1.0-SETUP.f_T-SETUP.f_C;

% topography
SETUP.e_R = 0.4;%0.2..0.4
SETUP.e_T = SETUP.e_R-0.1;

% hydraulic conductivities
SETUP.K_Reservoir = 5e-5;
SETUP.K=1e-5;

% stratigraphy
OL1 = [ 0.85, 0.00, 0.15, 1, 0.85 ] ;
OL2 = [ 0.75, 0.10, 0.15, 1, 0.75 ] ;
ML1 = [ 0.65, 0.30, 0.05, 2, 0.65 ] ;
EL1 = [ 0.65, 0.30, 0.05, 1, 0.55 ] ;
EL2 = [ 0.75, 0.20, 0.05, 1, 0.55 ] ;
EL3 = [ 0.85, 0.10, 0.05, 1, 0.55 ] ;
BL1 = [ 0.10, 0.90, 0.00, 1, 0.10 ] ;

stratigraphyMap = containers.Map( {'CENTER', 'RIM', 'TROUGH'}, ...
    { [ 0.00, OL1;...
        0.20, ML1;...
        0.80, EL1;...
        9.00, BL1 ],...
      [ 0.00, OL2;...
        0.10, ML1;...
        0.70, EL1;...
        0.90, EL2;...
        9.00+SETUP.e_R, BL1 ],...
      [ 0.00, OL1;
        0.20, ML1;
        0.40, EL2;
        0.60, EL3;
        9.00+SETUP.e_T, BL1 ] } );
    
SETUP.stratigraphy = { stratigraphyMap('CENTER'), ...
    stratigraphyMap('RIM'), ...
    stratigraphyMap('TROUGH') };

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
SETUP.runName = sprintf( [ 'SCENARIO_' SETUP.scenario '_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_xice%d_xH%d_xW%d_xS%d_%s_eRes%0.2f_snowDens%d_maxSnow%0.2f' ], ...
     SETUP.xice, SETUP.xH, SETUP.xW, SETUP.xS, SETUP.boundaryCondition_T, SETUP.e_Reservoir, SETUP.snowDens, SETUP.relMaxSnow ) ;
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
