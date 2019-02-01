% author: jnitzbon
% description: this script is used to specify individual setups for
% parallel runs of CryoGrid3. it is called by the script "submit_SLURM.sh"
% which is used to submit MATLAB jobs to the SLURM queue
clear all;
close all;

%% generate SETUP struct

SETUP = {};

SETUP.startFromRun = 0;     % set this to 0 for first run and to 1 if run needs to be continued from saved state

if SETUP.startFromRun


    SETUP.startFromRunDir='/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_mpi_polygon/runs';
    SETUP.startFromRunName='SCENARIO_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_geometry1_eRes-0.20_snowDens250_maxSnow0.40_CONTINUED_CONTINUED';
    SETUP.startFromRunYear=2083;

    % output directory
    SETUP.saveDir = SETUP.startFromRunDir;
    SETUP.runName = [ SETUP.startFromRunName '_CONTINUED' ];

    SETUP.numRealizations = 3; % needed to start parallel pool

else

    % parameters
    SETUP.numRealizations = 3;
    SETUP.syncTimestep=6./24;
    SETUP.startDate = datenum( 1999, 10, 1 );
    SETUP.endDate = datenum( 2099, 12, 31);
    SETUP.xH=1;
    SETUP.xW=1;
    SETUP.xS=1;
    SETUP.xE=1;
    SETUP.xice=1;
    
    SETUP.polygon_geometry = 2; % 1: hexagonal , 2: circular, cross-section

    SETUP.fieldCapacity = 0.50; % 0.40
    SETUP.relMaxSnow = 0.40; % 1.0
    SETUP.snowDens = 250;%200..250
    SETUP.boundaryCondition_T = 'DarcyReservoir';
    SETUP.e_Reservoir = -2.0;%-1.0;

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
    EL1 = [ 0.65, 0.30, 0.05, 1, 0.55 ] ;   % this is supposed to reflect only segregated and pore ice, no wedge ice
    EL2 = [ 0.75, 0.20, 0.05, 1, 0.55 ] ;
    EL3 = [ 0.85, 0.10, 0.05, 1, 0.55 ] ;
    EL4 = [ 0.95, 0.00, 0.05, 1, 0.55 ] ;
    BL1 = [ 0.10, 0.90, 0.00, 1, 0.10 ] ;   % bedrock layer
    
    polygonType = 'epigenetic'; %'epigenetic'

    stratigraphyMap = containers.Map( {'CENTERsyngenetic', 'RIMsyngenetic', 'TROUGHsyngenetic', 'CENTERepigenetic', 'RIMepigenetic', 'TROUGHepigenetic'}, ...
        { [ 0.00, OL1;...
            0.20, ML1;...
            0.90, EL1;...
            9.00, BL1 ],...
          [ 0.00, OL2;...
            0.10, ML1;...
            0.80, EL2;...
            9.00+SETUP.e_R, BL1 ],...
          [ 0.00, OL1;
            0.20, ML1;
            0.50, EL2;
            0.70, EL4;
            9.00+SETUP.e_T, BL1 ],...
          [ 0.00, OL1;...
            0.20, ML1;...
            0.90, EL1;...
            9.00, BL1 ],...
          [ 0.00, OL2;...
            0.10, ML1;...
            0.80, EL1;...
            9.00+SETUP.e_R, BL1 ],...
          [ 0.00, OL1;
            0.20, ML1;
            0.50, EL2;
            0.70, EL4;  % decreasing excess ice with depth reflecting narrowing ice-wedge
            1.20, EL3;
            1.70, EL2;
            2.20, EL1;  
            9.00+SETUP.e_T, BL1 ] } );

    SETUP.stratigraphy = { stratigraphyMap([ 'CENTER' polygonType ]), ...
        stratigraphyMap([ 'RIM' polygonType ]), ...
        stratigraphyMap([ 'TROUGH' polygonType ]) };

    SETUP.scenario='rcp85';

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
    %SETUP.runName = sprintf( [ 'SCENARIO_' SETUP.scenario '_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' ) '_xice%d_xH%d_xW%d_xS%d_%s_eRes%0.2f_snowDens%d_maxSnow%0.2f' ], ...
     %    SETUP.xice, SETUP.xH, SETUP.xW, SETUP.xS, SETUP.boundaryCondition_T, SETUP.e_Reservoir, SETUP.snowDens, SETUP.relMaxSnow ) ;
    SETUP.runName = sprintf( [ 'SCENARIO_' datestr( SETUP.startDate, 'yyyymm' ) '-' datestr(SETUP.endDate, 'yyyymm' )  '_' SETUP.scenario '_xice%d_xE%d_xH%d_xW%d_xS%d_%s_geometry%d_eRes%0.2f_snowDens%d_maxSnow%0.2f' ], ...
         SETUP.xice, SETUP.xE, SETUP.xH, SETUP.xW, SETUP.xS, polygonType, SETUP.polygon_geometry, SETUP.e_Reservoir, SETUP.snowDens, SETUP.relMaxSnow ) ; 
    [~, SETUP.git_commit_hash] = system('git rev-parse HEAD');

end


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
