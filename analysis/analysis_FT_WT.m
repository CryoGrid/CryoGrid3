function analysis_FT_WT(numTiles,year,LoadData,LoadDataRef)  % indicate number of tiles, year, load data off/on  e.g. p(2)lot_FT_WT(5,2000,1) 
% plot water and frost tables for all tiles for different months  
% set loadDataRef to 1, if reference file is in different Experiment directory (which has to be indicated...)

% do first plot_FT_WT.m to get zsoil and idx10m, idx20
% make sure that no variables are in workspace before script is run! (whole workspace is saved)
%ccc loadData switch needed? prob not...

if(LoadData==1)
%    ExpLoc='Prudhoe'; 
%    % ExpSet='exW2mm_snow4_SinS_xice'
%    ExpSetRef='exW2mm_snow4_SinS_xice'
%    ExpSet='exW2mm_snow4'
%     %ExpSetRef='exW2mm_snow4'
%    Scen='RCP85';
%    Runs='Runs_ERL_submission/';

    ExpLoc='Norilsk'
    ExpSet='ms190'
    ExpSetRef='ms190'
    Scen='ERA5';
    Runs='Runs_Norilsk/';
    Out=['N:/permarisk/CryoGrid3/',Runs];  
    %Out='M:/Norilsk/OUT/'
    OutDir=[Out,ExpSet,'/',num2str(numTiles),'tiles/']   
% single referene run (single tundra tile)    
    OutDirRef=[Out,ExpSetRef,'/1tiles/'];   % reference run for single tundra tile
    paraFileRef=[OutDirRef,ExpLoc,'_',Scen,'_',ExpSetRef,'_T1_settings']; 
    outFileRef= [OutDirRef,ExpLoc,'_',Scen,'_',ExpSetRef,'_out',num2str(year),'_T1']
    paraFile=[OutDir,ExpLoc,'_',Scen,'_',ExpSet,'_T1_settings']; outFile= [OutDir,ExpLoc,'_',Scen,'_',ExpSet,'_out',num2str(year),'_T1']
    if(LoadDataRef)
        load(paraFileRef); load(outFileRef)   
    else
        load(paraFile); load(outFile)           
    end
    Tsoil_tiles(:,:,numTiles+1) = OUT.cryoGrid3(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:); % last column: reference tile!
    LWC_tiles(:,:,numTiles+1)   = OUT.liquidWater(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:);
    WC_tiles(:,:,numTiles+1)    = OUT.water(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:);
    Qlat_tiles(:,:,numTiles+1)  = OUT.EB.Q_lateral(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:)/ (PARA.technical.syncTimeStep.*24.*3600); 
        % tile runs
    load(paraFile); load(outFile)   
    %zsoil_tiles = zeros(length(GRID.soil.soilGrid)-1,numTiles);
    ts = OUT.timestamp(); 
%    Tsoil_tiles = zeros(length(GRID.soil.soilGrid)-1,length(ts),numTiles+1); 
    Tsoil_10m_tiles = zeros(length(ts),numTiles); Tsoil_20m_tiles = zeros(length(ts),numTiles); Tsoil_z1_tiles = zeros(length(ts),numTiles); Tsoil_z2_tiles = zeros(length(ts),numTiles); 
    FT_tiles = zeros(length(ts),numTiles); WT_tiles = zeros(length(ts),numTiles);  Soil_surf_tiles = zeros(length(ts),numTiles); 
    Qnet_tiles = zeros(length(ts),numTiles); %Qlat_tiles = zeros(length(GRID.general.cT_grid),length(ts),numTiles);
    Qlat_intZ_tiles = zeros(length(ts),numTiles);
    %Qlat_z_tiles = zeros(length(GRID.general.cT_grid),length(ts),numTiles);  Qlat_intz_tiles = zeros(length(ts),numTiles);
    %dr_latWater = zeros(1,numTiles); 
    Wlat_tiles = zeros(length(ts),numTiles);
    %OUT.lateral.water_fluxes = zeros(numTiles,numTiles,length(ts));
    rain_tiles = zeros(length(ts),numTiles); snow_tiles = zeros(length(ts),numTiles);  evap_tiles = zeros(length(ts),numTiles); sublim_tiles = zeros(length(ts),numTiles);
    snowHeight_tiles = zeros(length(ts),numTiles); 
    
    zsoil_tiles = PARA.ensemble.initial_altitude - GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
    for n=1:numTiles
        paraFile=[OutDir,ExpLoc,'_',Scen,'_',ExpSet,'_T',num2str(n),'_settings'];  
        load(paraFile)
        outFile= [OutDir,ExpLoc,'_',Scen,'_',ExpSet,'_out',num2str(year),'_T',num2str(n)] 
        load(outFile)   
    %     T = OUT.cryoGrid3; eval(['T_',num2str(n),'=T;']);
        if(n==1); ts = OUT.timestamp(); time = datetime(ts,'ConvertFrom','datenum'); end
        %zsoil_tiles(:,n) = PARA.ensemble.initial_altitude(n) - GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
        [~,idx10] = min(abs(zsoil_tiles(:,n)+10)); [~,idx20] = min(abs(zsoil_tiles(:,n)+20));
         Tsoil_tiles(:,:,n) = OUT.cryoGrid3(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:);
         Tsoil_10m_tiles = squeeze(Tsoil_tiles(idx10,:,:)); Tsoil_20m_tiles = squeeze(Tsoil_tiles(idx20,:,:)); 
         Tsoil_z1_tiles = squeeze(Tsoil_tiles(1,:,:)); Tsoil_z2_tiles = squeeze(Tsoil_tiles(2,:,:)); 
         LWC_tiles(:,:,n) = OUT.liquidWater(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:);
         WC_tiles(:,:,n) = OUT.water(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:);
         Qlat_tiles(:,:,n) = OUT.EB.Q_lateral(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:) / (PARA.technical.syncTimeStep.*24.*3600);  

        % frost and water table    
        frost_table = OUT.location.infiltration_altitude(); water_table = OUT.location.water_table_altitude(); soil_level = OUT.location.soil_altitude();
        %eval(['FT_',num2str(n),'=frost_table;']); eval(['WT_',num2str(n),'=water_table;']);
        FT_tiles(:,n) = frost_table; WT_tiles(:,n) = water_table; Soil_surf_tiles(:,n) = soil_level;
    % net fluxes SEB
        Qnet_tiles(:,n) = OUT.EB.Qnet;  % W/m2       
    % lateral fluxes
       % if strcmp(ExpSet,'H1W1')
           % latHeat_tiles(:,n) = nansum(OUT.EB.Q_lateral);  % J/m2 per sync interval % OUT.EB.Q_lateral was calculated incorrectly before...!  OUT.EB.Q_lateral (z,t)
           %PARA.technical.syncTimeStep=PARA.technical.outputTimestep;
            Qlat_tiles(:,:,n) = OUT.EB.Q_lateral(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:) / (PARA.technical.syncTimeStep.*24.*3600);  
            Qlat_intZ_tiles(:,n) = squeeze(trapz(GRID.general.cT_grid,OUT.EB.Q_lateral)) / (PARA.technical.syncTimeStep.*24.*3600);  % in W/m2 - sum over depth
            %dE_tot_tiles(:,n) = nansum(OUT.lateral.dE_tot,2);  
            %%%Qlat_z_tiles(:,:,n) = squeeze(nansum(OUT.lateral.dE_cell,2)) /(PARA.technical.syncTimeStep *24.*3600);  % lat heat flux to each worker in W/m2  (z,tiles,time)  OUT.lateral.dE_cell (z,tiles,time)   sum over workers
            %%%Qlat_intz_tiles(:,n)= squeeze(trapz(GRID.general.cT_grid,Qlat_z_tiles(:,:,n)));
            Wlat_tiles(:,n) = OUT.WB.dr_lateralWater;
            %OUT.lateral.water_fluxes
       % end
    % water balance
        rain_tiles(:,n) = OUT.WB.dp_rain; snow_tiles(:,n) = OUT.WB.dp_snow; 
        evap_tiles(:,n) = OUT.WB.de; sublim_tiles(:,n) = OUT.WB.ds;
    % snow
        snowHeight_tiles(:,n) = OUT.snow.topPosition;
    end
    snowHeight_tiles(isnan(snowHeight_tiles))=0.; %replace NANs by 0

%    sum_WB_tiles = rain_tiles + snow_tiles + evap_tiles + sublim_tiles;
    sum_WB_tiles = rain_tiles + snowHeight_tiles/5. + evap_tiles + sublim_tiles; % todotodo use snow density instead of factor 5! use snow height instead of snowfall (accounting for limiting snow height to snowMax!) factor 1/5 accounts for densitiy diff between snow and water
    % find indeces of months (April 1st, May 1st, June 1st, Jul 1st, Aug 1st,...
    checkdate=datevec(ts);
    ind_Apr1=find(checkdate(:,2)==4 & checkdate(:,3)==1 & checkdate(:,4)==0 ); ind_May1=find(checkdate(:,2)==5 & checkdate(:,3)==1 & checkdate(:,4)==0 ); ind_Jun1=find(checkdate(:,2)==6 & checkdate(:,3)==1 & checkdate(:,4)==0 );
    ind_Jul1=find(checkdate(:,2)==7 & checkdate(:,3)==1 & checkdate(:,4)==0 ); ind_Aug1=find(checkdate(:,2)==8 & checkdate(:,3)==1 & checkdate(:,4)==0 ); ind_Sep1=find(checkdate(:,2)==9 & checkdate(:,3)==1 & checkdate(:,4)==0 );

    %% Define soil layers. #1 Vertice-coordinates #2 Polygon side numbers #3 Color (RGB 255)
    % bedrock = [0 -11; 0 -10; 25 -10; 25 -11];  bedrockFaces = [1 2 3 4];  bedrockColors = [197 197 197]./255;
    mineralsoil = [0 -10; 0 -1.5; 10 -1.5; 10 -0.3; 50 -0.3; 50 -10]; mineralsoilFaces = [1 2 3 4 5 6]; mineralsoilColors = [216 201 180]./255;
    gravel = [0 -1.5; 0 2.5; 5 2.5; 10 0; 10 -1.5]; gravelFaces = [1 2 3 4 5]; gravelColors = [98 97 114]./255;
    peat = [10 -0.3; 10 0; 50 0; 50 -0.3]; peatFaces = [1 2 3 4]; peatColors = [128 122 83]./255;
%    snow = [5 2.5; 15 0.4; 25 0.4; 25 0; 10 0]; snowFaces = [1 2 3 4 5]; snowColors = [200 237 248]./255;
    snow = [5 2.5; 15 0.3; 50 0.3; 50 0; 10 0]; snowFaces = [1 2 3 4 5]; snowColors = [200 237 248]./255;   % set to 0.4 cccc

    % PARA.ensemble.TileWidth(end)=PARA.ensemble.TileWidth(end-1);
    % disp('attention: last tile width set to previous!')
    distanceRC=cumsum(PARA.ensemble.TileWidth)-0.5*PARA.ensemble.TileWidth; % distance to Road Centre in meters
    tileBoarders=cumsum(PARA.ensemble.TileWidth); % distance of upper tile boarders to Road Centre in metres
    % %x_GRprofile = [0 tileBoarders]; 
    % %x_GRprofile = 0:5:25; y_GRprofile = [PARA.IS.EHag PARA.IS.EHag 0 0 0 0];
    % x_GRprofile = 0:5:25; y_GRprofile = [PARA.IS.EHag PARA.IS.EHag PARA.IS.EHag/2.  0 0 0];    
    
    save(['Data/',Runs,'/Data_',num2str(numTiles),'tiles_',num2str(year),'_',ExpSet])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
z1=1; z2=110;
%zzz=squeeze(T_tiles(:,ind_Aug1,:));
imagesc(1:5,zsoil_tiles(z1:z2,1),squeeze(Tsoil_tiles(z1:z2,ind_Aug1,:)))     
set(gca,'YDir','normal');  title(['1st of August ',num2str(year),' ',ExpLoc,ExpSet],'Interpreter','none'); 
colorbar

figure(11)
z1=1; z2=110;
%zzz=squeeze(T_tiles(:,ind_Aug1,:));
imagesc(1:5,zsoil_tiles(z1:z2,1),squeeze(LWC_tiles(z1:z2,ind_Aug1,:)))     
set(gca,'YDir','normal');  title(['1st of August ',num2str(year),' ',ExpLoc,ExpSet],'Interpreter','none'); 
colorbar

figure(2)
cmap=winter(numTiles);
    subplot(2,1,1)
hold on
for n=1:numTiles
    h(n)=plot(squeeze(Tsoil_tiles(z1:z2,ind_Aug1,n)),zsoil_tiles(z1:z2,n),'Color',cmap(n,:));
    title(['soil temperatures at year ',num2str(year)])
end
hold off; grid on; 
    subplot(2,1,2)
hold on
for n=1:numTiles
    h1=plot(time,Tsoil_10m_tiles(:,n),'Color',cmap(n,:));
    h2=plot(time,Tsoil_20m_tiles(:,n),'--','Color',cmap(n,:));
end
grid on; legend([h1(1) h2(1)],'Tsoil 10m','Tsoil 20m')

figure(3) % plot frost and water tables for all tiles
    subplot(2,2,1)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Silty soil');
hold on
%plot(distanceRC,FT_tiles(ind_Jun1,:),'+-r',distanceRC,WT_tiles(ind_Jun1,:),'+-b',tileBoarders,Soil_surf_tiles(ind_Jun1,:),'+--k')
%plot(distanceRC,FT_tiles(ind_Jun1,:),'+-r',distanceRC,WT_tiles(ind_Jun1,:),'+-b',distanceRC,Soil_surf_tiles(ind_Jun1,:),'+--k',x_GRprofile,y_GRprofile,'m')
plot(distanceRC,FT_tiles(ind_Jun1,:),'+-r','DisplayName','Frost table')
plot(distanceRC,WT_tiles(ind_Jun1,:),'+-b','DisplayName','Water table') %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
%plot(distanceRC,snowHeight_tiles(ind_Jun1,:),'+-c','DisplayName','Snow height')
plot(distanceRC,nanmax(snowHeight_tiles+Soil_surf_tiles(1,:)),'+-c','DisplayName','Max Snow height') % maximum snow height in year

% hold on; plot(distanceRC(end-1),FT_tiles(ind_Jun1,end),'*r',distanceRC(end-1),WT_tiles(ind_Jun1,end),'*b',tileBoarders(end-1),Soil_surf_tiles(ind_Jun1,end),'*k'); hold off
%grid on; axis([0 numTiles -1 2.5]); %%grid on; axis([0 distanceRC(end-1) -1 2.5])
grid on; axis([0 tileBoarders(end) -8 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of June ',num2str(year),' ',ExpLoc,ExpSet],'Interpreter','none'); 
legend('show', 'Location', 'southeast'); % legend('frost table','water table','soil surface')
    subplot(2,2,2)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Silty soil');
hold on
plot(distanceRC,FT_tiles(ind_Jul1,:),'+-r','DisplayName','Frost table')
plot(distanceRC,WT_tiles(ind_Jul1,:),'+-b','DisplayName','Water table') 
grid on; axis([0 tileBoarders(end) -8 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of July ',num2str(year),' ',ExpLoc,ExpSet],'Interpreter','none'); 
legend('show', 'Location', 'southeast'); 
    subplot(2,2,3)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Silty soil');
hold on
plot(distanceRC,FT_tiles(ind_Aug1,:),'+-r','DisplayName','Frost table')
plot(distanceRC,WT_tiles(ind_Aug1,:),'+-b','DisplayName','Water table') 
grid on; axis([0 tileBoarders(end) -8 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of Aug ',num2str(year),' ',ExpLoc,ExpSet],'Interpreter','none'); 
legend('show', 'Location', 'southeast'); 
    subplot(2,2,4)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Silty soil');
hold on
plot(distanceRC,FT_tiles(ind_Sep1,:),'+-r','DisplayName','Frost table')
plot(distanceRC,WT_tiles(ind_Sep1,:),'+-b','DisplayName','Water table') 
grid on; axis([0 tileBoarders(end) -8 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of Sep ',num2str(year),' ',ExpLoc,ExpSet],'Interpreter','none'); 
legend('show', 'Location', 'southeast'); 



figure(4) % z-integrated lateral heat fluxes
    subplot(1,2,1)
hold on
h1=plot(time,Qlat_intZ_tiles,time,sum(Qlat_intZ_tiles,2),'m');
grid on; legend(h1,'T1','T2','T3','T4','T5','sum'); title('z integrated lateral heat fluxes') %todotodo: sum is not zero, this would be only the case if all tiles on same surface level...
hold off
    subplot(1,2,2)
plot(time,Qnet_tiles,'c')
hold on
h1=plot(time,Qlat_intZ_tiles,time,sum(Qlat_intZ_tiles,2),'m');
grid on; legend(h1,'T1','T2','T3','T4','T5','sum'); title('Qnet and z integrated lateral heat fluxes')
hold off

% figure(5) % vertically resolved lateral heat fluxes
%     subplot(2,2,1)
% plot(1:numTiles,squeeze(Qlat_tiles(:,ind_Jun1,:)),'+-r')
% grid on
% xlabel('Tile number'); ylabel('lat heat [W/m2]'); title(['1st of June',num2str(year)])
%     subplot(2,2,2)
% plot(1:numTiles,Qlat_tiles(ind_Jul1,:),'+-r')
% grid on
% xlabel('Tile number'); ylabel('lat heat [W/m2]'); title(['1st of July',num2str(year)])
%     subplot(2,2,3)
% plot(1:numTiles,Qlat_tiles(ind_Aug1,:),'+-r')
% grid on
% xlabel('Tile number'); ylabel('lat heat [W/m2]'); title(['1st of August',num2str(year)])
%     subplot(2,2,4)
% plot(1:numTiles,Qlat_tiles(ind_Sep1,:),'+-r')
% grid on
% xlabel('Tile number'); ylabel('lat heat [W/m2]'); title(['1st of September',num2str(year)])
    

figure(6) % Qnet
    subplot(2,2,1)
plot(1:numTiles,Qnet_tiles(ind_Jun1,:),'+-r')
grid on
xlabel('Tile number'); ylabel('Qnet [W/m2]'); title(['1st of June',num2str(year)])
    subplot(2,2,2)
plot(1:numTiles,Qnet_tiles(ind_Jul1,:),'+-r')
grid on
xlabel('Tile number'); ylabel('Qnet [W/m2]'); title(['1st of July',num2str(year)])
    subplot(2,2,3)
plot(1:numTiles,Qnet_tiles(ind_Aug1,:),'+-r')
grid on
xlabel('Tile number'); ylabel('Qnet [W/m2]'); title(['1st of August',num2str(year)])
    subplot(2,2,4)
plot(1:numTiles,Qnet_tiles(ind_Sep1,:),'+-r')
grid on
xlabel('Tile number'); ylabel('Qnet [W/m2]'); title(['1st of September',num2str(year)])


figure(7) % lateral water fluxes for all workers
    subplot(2,2,1)
ti=ind_Jun1;
plot(1:numTiles,Wlat_tiles(ti,:),'+-r')
grid on; title(['lateral water Jun1   sum=',num2str(nansum(Wlat_tiles(ti,:),2))])
    subplot(2,2,2)
ti=ind_Jul1;
plot(1:numTiles,Wlat_tiles(ti,:),'+-r')
grid on; title(['lateral water Jul1   sum=',num2str(nansum(Wlat_tiles(ti,:),2))])
    subplot(2,2,3)
ti=ind_Aug1;
plot(1:numTiles,Wlat_tiles(ti,:),'+-r')
grid on; title(['lateral water Aug1   sum=',num2str(nansum(Wlat_tiles(ti,:),2))])
    subplot(2,2,4)
ti=ind_Sep1;
plot(1:numTiles,Wlat_tiles(ti,:),'+-r')
grid on; title(['lateral water Sep1    sum=',num2str(nansum(Wlat_tiles(ti,:),2))])

    figure(8) % water balance
plot(time,cumsum(rain_tiles(:,numTiles)),'b',time,cumsum(snowHeight_tiles(:,numTiles)/5.),'c',time,cumsum(snow_tiles(:,numTiles)),'c--',time,-cumsum(evap_tiles(:,numTiles)),'r',time,-cumsum(sublim_tiles(:,numTiles)),'g',time,cumsum(sum_WB_tiles(:,numTiles)),'m') %todotodo  /5
grid on; title(['water balance Tile ',num2str(numTiles)]); legend('rain','snowHeight/5','snow','-evap','-sublim','sum');

    figure(9) % snow height
plot(time,snowHeight_tiles)
title('snow height'); grid on; legend('SU1','SU2','SU3','SU4a','SU4b')

end