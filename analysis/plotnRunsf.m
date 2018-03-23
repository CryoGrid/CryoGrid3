function Coupled = plotnRunsf( run_name )
% Function used to plot results from Leo's version (March 2017)

%% Adjust names
while  isnan(str2double(run_name(end)))==0;
    run_name(end)=[];
end

close all
clearvars -except run_name

%% Load data and display degree day statistics
addpath \\lagringshotell\geofag\projects\coup\Leo\toolbox
Coupled(100).OUT='Gerard';
nbreal=1;
while 7==exist([run_name num2str(nbreal)],'dir') == 1;
    filename=[run_name num2str(nbreal)];
    Names=dir(['.\' filename]);
    dirFlags=[Names.isdir];
    Names=Names(logical(1-dirFlags));
    Names={Names.name}';
    i_final=1;
    i_out=1;
    for i=1:length(Names);
        if  isempty(strfind(Names{i},'final'))==0;
            Names_final{i_final,1}=Names{i};
            i_final=i_final+1;
        elseif isempty(strfind(Names{i},'output'))==0;
            Names_out{i_out,1}=Names{i};
            i_out=i_out+1;
        else
            Names_settings=Names{i};
        end
    end
    n_stack=5;% Number of years stacked ---------------------------<<<<
    if n_stack>length(Names_out);
        n_stack=length(Names_out);
    end
    for i=1:n_stack;
        name_out = fullfile(filename,Names_out{end-n_stack+i});
        name_fin = fullfile(filename,Names_final{end-n_stack+i});
        load(name_out)
        load(name_fin)
        OUTf(i)=OUT;
        FINALf(i)=FINAL;
    end
    Coupled(nbreal).OUT=OUTf;
    Coupled(nbreal).FINAL=FINALf;
    name_set = fullfile(filename,Names_settings);
    load(name_set)
    Coupled(nbreal).GRID=GRID;
    Coupled(nbreal).PARA=PARA;
    nbreal=nbreal+1;
    DDf_1yr( '.',filename,1 );
end
Coupled(nbreal:end)=[];
nbreal=nbreal-1;

clear GRID PARA OUT filename name_out name_set run_number dirFlags Names ans

%% Plot all data in 1 figure


% Merge data to plot
for real=1:nbreal;
    Coupled(real).PLOT.timestamp=[];
    Coupled(real).PLOT.cryoGrid3=[];
    Coupled(real).PLOT.water_table_altitude=[];
    Coupled(real).PLOT.active_layer_depth_altitude=[];
    Coupled(real).PLOT.dr_lateral=[];
    Coupled(real).PLOT.dr_DarcyReservoir=[];
    Coupled(real).PLOT.dr_lateralExcess=[];
    for year=1:length(Coupled(real).OUT)
        Coupled(real).PLOT.timestamp=[Coupled(real).PLOT.timestamp ; Coupled(real).OUT(year).timestamp];
        Coupled(real).PLOT.cryoGrid3=[Coupled(real).PLOT.cryoGrid3 Coupled(real).OUT(year).cryoGrid3];
        Coupled(real).PLOT.water_table_altitude=[Coupled(real).PLOT.water_table_altitude ; Coupled(real).OUT(year).location.water_table_altitude];
        Coupled(real).PLOT.active_layer_depth_altitude=[Coupled(real).PLOT.active_layer_depth_altitude ; Coupled(real).OUT(year).location.active_layer_depth_altitude];
        Coupled(real).PLOT.dr_lateral=[Coupled(real).PLOT.dr_lateral ; Coupled(real).OUT(year).WB.dr_lateral];
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
    plot(Coupled(i).PLOT.timestamp,Coupled(i).PLOT.active_layer_depth_altitude,'k')
    ylabel('Elevation (masl)')
    hold off
end

% Worker fluxes and boundary condition fluxes
for i=1:nbreal;
    subplot(3,nbreal,i+2*nbreal)
    % nonnul=Coupled(i).PLOT.dr_lateral~=0;
    plot(Coupled(i).PLOT.timestamp,cumsum(Coupled(i).PLOT.dr_lateral),'.')
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