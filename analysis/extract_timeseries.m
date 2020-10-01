% extract variables for timeseries plot for indicated tile number
clear all

numTiles=5
startyear=1990; endyear=2010;
years=startyear:endyear;
iii=0; % counter
%zLminTpos=zeros(length(years),numTiles); zLminLWC=zeros(length(years),numTiles); % depth z of Thaw Depth maximum (based on T>0 criterium)

SetLoc='Prudhoe'; 
ExpSet='exW2mm_snow4_SinS_xice'
%SetLoc='HapVal'
%ExpSet='exW2mm_snow4_new'
Scen='RCP85';
Out='N:/permarisk/CryoGrid3/Runs_ERL_submission/';    
%Out = '/data/permarisk/CryoGrid3/Runs_ERL_submission/'; 
OutDir=[Out,ExpSet,'/',num2str(numTiles),'tiles/']; %OutDir=[Out,ExpSet,'/',num2str(numTiles),'tiles/'];

paraFile=[OutDir,SetLoc,'_',Scen,'_',ExpSet,'_T1_settings']; 
outFile=[OutDir,SetLoc,'_',Scen,'_',ExpSet,'_out',num2str(startyear),'_T1']
load(paraFile); load(outFile) 
forcingFile=['../forcing/',PARA.forcing.filename];
load(forcingFile);

startyearF=find(ismember(FORCING.data.t_span,datenum(startyear,1,1))); endyearF=find(ismember(FORCING.data.t_span,datenum(endyear,1,1))); %find indices in forcing timevector to cover same spanned by startyear to endyear
z1=GRID.soil.cT_domain_ub; [~,z2]=min(abs(GRID.general.cT_grid-50)); % constrain soil domain to 50m
zsoil_tiles = PARA.ensemble.initial_altitude - GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
        
%% forcing
timeF=FORCING.data.time(startyearF:endyearF); TairF=FORCING.data.Tair(startyearF:endyearF); rainF=FORCING.data.rainfall(startyearF:endyearF); snowF=FORCING.data.snowfall(startyearF:endyearF);

%%
for n=1:numTiles  % loop over tiles
    [~,idx5(n)] = min(abs(zsoil_tiles(:,n)+5)); [~,idx10(n)] = min(abs(zsoil_tiles(:,n)+10)); [~,idx20(n)] = min(abs(zsoil_tiles(:,n)+20));
    for i=1:length(years) % year loop  
        outFile = [OutDir,SetLoc,'_',Scen,'_',ExpSet,'_out',num2str(years(i)),'_T',num2str(n)]  
        load(outFile)   
        if(iii==0)  % first year only      
            if(n==1); ts = OUT.timestamp'; end % outputtime vector is the same for all tiles
            Tsoil = OUT.cryoGrid3(z1:z2,:); eval(['Tsoil_T',num2str(n),'=Tsoil;']);
            LWC = OUT.liquidWater(z1:z2,:); eval(['LWC_T',num2str(n),'=LWC;']);
            WC = OUT.water(z1:z2,:); eval(['WC_T',num2str(n),'=WC;']);            
            %%% T_tiles(:,1:length(ts),n) = OUT.cryoGrid3; eval('T = OUT.cryoGrid3;');
            [zL1Tpos_ub,zL1Tpos_lb,zL2Tpos_ub,zL2Tpos_lb,zL3Tpos_ub,zL3Tpos_lb] = find_thawedLayers(zsoil_tiles(:,n),Tsoil,0); % thawed layer based on T
            [zL1LWC_ub,zL1LWC_lb,zL2LWC_ub,zL2LWC_lb,zL3LWC_ub,zL3LWC_lb] = find_thawedLayers(zsoil_tiles(:,n),LWC./WC,1);     % thawed layers based on LWC
            zLminTpos(i,n)=nanmin([zL1Tpos_lb,zL2Tpos_lb,zL3Tpos_lb]);
            zLminLWC(i,n)=nanmin([zL1LWC_lb,zL2LWC_lb,zL3LWC_lb]);            
            eval(['zL1Tpos_ub_T',num2str(n),'=zL1Tpos_ub;']); eval(['zL1Tpos_lb_T',num2str(n),'=zL1Tpos_lb;']);
            eval(['zL2Tpos_ub_T',num2str(n),'=zL2Tpos_ub;']); eval(['zL2Tpos_lb_T',num2str(n),'=zL2Tpos_lb;']);
            eval(['zL3Tpos_ub_T',num2str(n),'=zL3Tpos_ub;']); eval(['zL3Tpos_lb_T',num2str(n),'=zL3Tpos_lb;']);
            FT = OUT.location.infiltration_altitude; eval(['FT_T',num2str(n),'=FT;']); eval(['FTmin_T',num2str(n),'=min(FT);']); 
            WT = OUT.location.water_table_altitude;  eval(['WT_T',num2str(n),'=WT;']);
            SL = OUT.location.soil_altitude;       eval(['SL_T',num2str(n),'=SL;']);
            SnowHeight = OUT.snow.topPosition;     eval(['SnowHeight_T',num2str(n),'=SnowHeight;']);
            
            Qnet = OUT.SEB.QNET;                     eval(['Qnet_T',num2str(n),'=Qnet;']);
            Qh = OUT.SEB.QH;                         eval(['Qh_T',num2str(n),'=Qh;']);
            Qe = OUT.SEB.QE;                         eval(['Qe_T',num2str(n),'=Qe;']);
            
            QlatZ = squeeze(trapz(GRID.general.cT_grid,OUT.EB.Q_lateral))' / (PARA.technical.syncTimeStep.*24.*3600);  % in W/m2 - sum over depth
            eval(['QlatZ_T',num2str(n),'=QlatZ;']);
            
            rainOut = OUT.WB.dp_rain; eval(['rainOut_T',num2str(n),'=rainOut;']);
            snowOut = OUT.WB.dp_snow; eval(['snowOut_T',num2str(n),'=snowOut;']);
            iii=1;
        else
            if(n==1); eval('ts = [ts,OUT.timestamp''];'); end
            Tsoil = OUT.cryoGrid3(z1:z2,:); eval(['Tsoil_T',num2str(n),' = [Tsoil_T',num2str(n),' Tsoil];']);
            LWC = OUT.liquidWater(z1:z2,:); eval(['LWC_T',num2str(n),' = [LWC_T',num2str(n),' LWC];']);
            WC = OUT.water(z1:z2,:); eval(['WC_T',num2str(n),' = [WC_T',num2str(n),' WC];']);
            %eval(['Tsoil_T',num2str(n),' = [Tsoil_T',num2str(n),' OUT.cryoGrid3(z1:z2,:)];']);
            %%% zzz=squeeze(Tsoil_Tiles(:,:,n)); eval('Tsoil_Tiles(:,1:size(zzz,2),n) = cat(2,squeeze(Tsoil_Tiles(:,1:size(zzz_previous,2),n)),OUT.cryoGrid3);'); 
            [zL1Tpos_ub,zL1Tpos_lb,zL2Tpos_ub,zL2Tpos_lb,zL3Tpos_ub,zL3Tpos_lb] = find_thawedLayers(zsoil_tiles(:,n),Tsoil,0);                     
            [zL1LWC_ub,zL1LWC_lb,zL2LWC_ub,zL2LWC_lb,zL3LWC_ub,zL3LWC_lb] = find_thawedLayers(zsoil_tiles(:,n),LWC./WC,1);            
            zLminTpos(i,n)=min([zL1Tpos_lb,zL2Tpos_lb,zL3Tpos_lb]);
            zLminLWC(i,n)=min([zL1LWC_lb,zL2LWC_lb,zL3LWC_lb]);
            %eval(['zL1Tpos_ub_T',num2str(n),'= zL1Tpos_ub;']); eval(['zL1Tpos_lb_T',num2str(n),'=zL1Tpos_lb;']);      
            eval(['zL1Tpos_ub_T',num2str(n),' = [zL1Tpos_ub_T',num2str(n),' zL1Tpos_ub];']); eval(['zL1Tpos_lb_T',num2str(n),' = [zL1Tpos_lb_T',num2str(n),' zL1Tpos_lb];']);
            eval(['zL2Tpos_ub_T',num2str(n),' = [zL2Tpos_ub_T',num2str(n),' zL2Tpos_ub];']); eval(['zL2Tpos_lb_T',num2str(n),' = [zL2Tpos_lb_T',num2str(n),' zL2Tpos_lb];']);
            eval(['zL3Tpos_ub_T',num2str(n),' = [zL3Tpos_ub_T',num2str(n),' zL3Tpos_ub];']); eval(['zL3Tpos_lb_T',num2str(n),' = [zL3Tpos_lb_T',num2str(n),' zL3Tpos_lb];']);
            % same for zLxLWC!? ccc
            eval(['FT_T',num2str(n),' = [FT_T',num2str(n),'; OUT.location.infiltration_altitude];']); eval(['FTmin_T',num2str(n),' = [FTmin_T',num2str(n),'; min(OUT.location.infiltration_altitude)];']);
            eval(['WT_T',num2str(n),' = [WT_T',num2str(n),'; OUT.location.water_table_altitude];']);
            eval(['SL_T',num2str(n),' = [SL_T',num2str(n),'; OUT.location.soil_altitude];']);
            eval(['SnowHeight_T',num2str(n),' = [SnowHeight_T',num2str(n),'; OUT.snow.topPosition];']);
            eval(['Qnet_T',num2str(n),' = [Qnet_T',num2str(n),'; OUT.SEB.QNET];']); eval(['Qh_T',num2str(n),' = [Qh_T',num2str(n),'; OUT.SEB.QH];']);  eval(['Qe_T',num2str(n),' = [Qe_T',num2str(n),'; OUT.SEB.QE];']);            
            QlatZ = squeeze(trapz(GRID.general.cT_grid,OUT.EB.Q_lateral))' / (PARA.technical.syncTimeStep.*24.*3600);
            eval(['QlatZ_T',num2str(n),' = [QlatZ_T',num2str(n),'; QlatZ];']);
            eval(['rainOut_T',num2str(n),' = [rainOut_T',num2str(n),'; OUT.WB.dp_rain];']);
            eval(['snowOut_T',num2str(n),' = [snowOut_T',num2str(n),'; OUT.WB.dp_snow];']);
        end
    end
    %%%eval(['T_tiles(:,:,n)=T_T',num2str(n),';']) 
    %%%eval('Tsoil_tiles(:,:,n) = T_tiles(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:,n);')
    eval(['Tsoil_tiles(:,:,n)=Tsoil_T',num2str(n),';']) %cccc could be coded more efficiently without extracting _T fields for individual tiles...!?
    eval('Tsoil_5m_tiles = squeeze(Tsoil_tiles(idx5(n),:,:)); Tsoil_10m_tiles = squeeze(Tsoil_tiles(idx10(n),:,:)); Tsoil_20m_tiles = squeeze(Tsoil_tiles(idx20(n),:,:));')
    eval(['LWC_tiles(:,:,n)=LWC_T',num2str(n),';'])
    eval(['FT_tiles(:,n)=FT_T',num2str(n),';']);  eval(['WT_tiles(:,n)=WT_T',num2str(n),';']);
    eval(['SL_tiles(:,n)=SL_T',num2str(n),';']); 
    eval(['SnowHeight_tiles(:,n)=SnowHeight_T',num2str(n),';']); 
    eval(['Qnet_tiles(:,n)=Qnet_T',num2str(n),';']); eval(['Qh_tiles(:,n)=Qh_T',num2str(n),';']);  eval(['Qe_tiles(:,n)=Qe_T',num2str(n),';']); 
    eval(['QlatZ_tiles(:,n)=QlatZ_T',num2str(n),';']); 
    eval(['rainOut_tiles(:,n)=rainOut_T',num2str(n),';']);  eval(['WT_tiles(:,n)=WT_T',num2str(n),';']);
    eval(['snowOut_tiles(:,n)=snowOut_T',num2str(n),';']);  eval(['WT_tiles(:,n)=WT_T',num2str(n),';']);
    iii=0;
end
saveDir='Data/Runs_ERL_submission';
%saveDir='Test';

if ~exist(saveDir,'dir'); mkdir(saveDir); end 
% Tsoil matrix gets rather large for long timeseries!
%save([saveDir,'/Data',num2str(numTiles),'_T_',SetLoc,'_',Scen,'_',ExpSet],'ts*','zsoil*','idx*','Tsoil*','-v7.3'); % data can be stored >2GB, but compression makes things slow... 
save([saveDir,'/Data',num2str(numTiles),'_Tsoil_10m20m_',SetLoc,'_',Scen,'_',ExpSet],'ts*','Tsoil_5m*','Tsoil_10m*','Tsoil_20m*'); 
save([saveDir,'/Data',num2str(numTiles),'_FTWT_',SetLoc,'_',Scen,'_',ExpSet],'ts*','zsoil*','FT*','zL*','WT*','SL*','SnowHeight*');
save([saveDir,'/Data',num2str(numTiles),'_Q_',SetLoc,'_',Scen,'_',ExpSet],'ts*','zsoil*','Q*');
save([saveDir,'/Data',num2str(numTiles),'_F_',SetLoc,'_',Scen],'*rain*','*snow*','timeF','TairF')  
%save([saveDir,'/Data',num2str(numTiles),'_LWC_',SetLoc,'_',Scen,'_',ExpSet],'ts*','zsoil*','idx*','LWC*','-v7.3');
