function Coupled = plotnRunsf( run_name )
% Function used to plot results from Leo's version (March 2017)

assert(7==exist(run_name,'dir'),'The folder does not exist chico')

close all
clearvars -except run_name

%% Load data and display degree day statistics

% Retrieve file content
addpath \\lagringshotell\geofag\projects\coup\Leo\toolbox
[ name_str ] = runFiles( ['.\' run_name] );
nbreal=length(name_str);
nbyear=length(name_str(1).output);
Coupled(nbreal).OUT='Gerard';

% Prepare stacking
n_stack=1;% Number of years stacked ---------------------------<<<<
if n_stack>nbyear;
    n_stack=nbyear;
end

for ireal=1:nbreal
    
    % Stack
    for i=1:n_stack;
        name_out = fullfile(run_name,name_str(ireal).output{end-n_stack+i});
        name_fin = fullfile(run_name,name_str(ireal).final{end-n_stack+i});
        load(name_out)
        load(name_fin)
        OUTf(i)=OUT;
        FINALf(i)=FINAL;
    end
    
    % Store
    Coupled(ireal).OUT=OUTf;
    Coupled(ireal).FINAL=FINALf;
    name_set = fullfile(run_name,name_str(ireal).settings);
    load(name_set)
    Coupled(ireal).GRID=GRID;
    Coupled(ireal).PARA=PARA;
    
    % Display DD values
    DDf_1yr_direct( run_name, name_str(ireal).output{end}, name_str(ireal).settings, 1 );
    
end

clear FORCING GRID PARA FINAL FINALf OUT OUTf i ireal name_out name_set name_fin ans

%% Plot all data in 1 figure


% Merge data to plot
for real=1:nbreal;
    Coupled(real).PLOT.timestamp=[];
    Coupled(real).PLOT.cryoGrid3=[];
    Coupled(real).PLOT.water_table_altitude=[];
    Coupled(real).PLOT.infiltration_altitude=[];
    Coupled(real).PLOT.dr_lateralWater=[];
    Coupled(real).PLOT.dr_DarcyReservoir=[];
    Coupled(real).PLOT.dr_lateralExcess=[];
    for year=1:length(Coupled(real).OUT)
        Coupled(real).PLOT.timestamp=[Coupled(real).PLOT.timestamp ; Coupled(real).OUT(year).timestamp];
        Coupled(real).PLOT.cryoGrid3=[Coupled(real).PLOT.cryoGrid3 Coupled(real).OUT(year).cryoGrid3];
        Coupled(real).PLOT.water_table_altitude=[Coupled(real).PLOT.water_table_altitude ; Coupled(real).OUT(year).location.water_table_altitude];
        Coupled(real).PLOT.infiltration_altitude=[Coupled(real).PLOT.infiltration_altitude ; Coupled(real).OUT(year).location.infiltration_altitude];
        Coupled(real).PLOT.dr_lateralWater=[Coupled(real).PLOT.dr_lateralWater ; Coupled(real).OUT(year).WB.dr_lateralWater];
        Coupled(real).PLOT.dr_DarcyReservoir=[Coupled(real).PLOT.dr_DarcyReservoir ; Coupled(real).OUT(year).WB.dr_DarcyReservoir];
        Coupled(real).PLOT.dr_lateralExcess=[Coupled(real).PLOT.dr_lateralExcess ; Coupled(real).OUT(year).WB.dr_lateralExcess];
    end
end

figure('Position', [35, 50, 1850, 1050]);
% Crygrid figs
for i=1:nbreal;
    subplot(3,nbreal,i)
    graphtitle=[run_name num2str(i)];
    tirets= graphtitle == '_';
    graphtitle(tirets)='-';
    [nl,nc]=size(Coupled(i).PLOT.cryoGrid3);
    i_starttime=1; % nc-2920+1;
    i_startz=10;
    i_endz=115;
    ElevVect=(Coupled(i).PARA.location.altitude-Coupled(i).GRID.general.K_grid(i_startz:i_endz));
    
    imagesc(Coupled(i).PLOT.timestamp(i_starttime:end),ElevVect,Coupled(i).PLOT.cryoGrid3(i_startz+16:i_endz+16,i_starttime:end))
    hold on
    set(gca,'YDir','normal')
    title(graphtitle)
    colorbar
    contour(Coupled(i).PLOT.timestamp(i_starttime:end),ElevVect,Coupled(i).PLOT.cryoGrid3(i_startz:i_endz,i_starttime:end),[0 0],'black')
    set(gca,'YDir','normal')
    datetick
    ylabel('Elevation (masl)')
    hold off
end

% Water table and bucket figure
for i=1:nbreal;
    subplot(3,nbreal,i+nbreal)
    plot(Coupled(i).PLOT.timestamp,Coupled(i).PLOT.water_table_altitude,'.')
    hold on
    datetick
    title('Water table and bottom bucket (masl)')
    plot(Coupled(i).PLOT.timestamp,Coupled(i).PLOT.infiltration_altitude,'k')
    ylabel('Elevation (masl)')
    hold off
end

% Worker fluxes and boundary condition fluxes
for i=1:nbreal;
    subplot(3,nbreal,i+2*nbreal)
    % nonnul=Coupled(i).PLOT.dr_lateralWater~=0;
    plot(Coupled(i).PLOT.timestamp,cumsum(Coupled(i).PLOT.dr_lateralWater),'.')
    hold on
    datetick
    title('Network and BC water (mm water change)')
    ylabel('cumulative lateral water changes (mm for worker)')
    if nbreal>1
        if strcmp(Coupled(i).PARA.ensemble.boundaryCondition(i).type, 'DarcyReservoir')==1;
            % nonnul=Coupled(i).PLOT.dr_DarcyReservoir~=0;
            plot(Coupled(i).PLOT.timestamp,cumsum(Coupled(i).PLOT.dr_DarcyReservoir))
        else
            % nonnul=Coupled(i).PLOT.dr_lateralExcess~=0;
            plot(Coupled(i).PLOT.timestamp,cumsum(Coupled(i).PLOT.dr_lateralExcess))
        end
    end
    hold off
end

end