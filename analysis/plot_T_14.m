% plot T timeseries for various depth

%filespec='MSW_xH1_xW0_xS0_infil1_xice0_rF1_sF1';
%filespec='MSW_xH1';
filespec='LF25_xH1'
yearspec='2010';

% NL
outputfile1 = ['../runs/WD0_',filespec,'/WD0_',filespec,'_',yearspec];
configfile1 = ['../runs/WD0_',filespec,'/WD0_',filespec,'_settings'];
% Lake
outputfile2 = ['../runs/WD5_',filespec,'/WD5_',filespec,'_',yearspec];
configfile2 = ['../runs/WD5_',filespec,'/WD5_',filespec,'_settings'];

% % NL
% outputfile1 = ['../runs/',filespec,'_i1/',filespec,'_i1_output',yearspec];
% configfile1 = ['../runs/',filespec,'_i1/',filespec,'_i1_settings'];
% %Lake
% outputfile2 = ['../runs/',filespec,'_i2/',filespec,'_i2_output',yearspec];
% configfile2 = ['../runs/',filespec,'_i2/',filespec,'_i2_settings'];


% outputfile1 = ['../runs/noLateral/runs/',filespec,'_i2/',filespec,'_i2_output',yearspec];
% configfile1 = ['../runs/noLateral/runs/',filespec,'_i2/',filespec,'_i2_settings'];
% % Lake
% outputfile2 = ['../runs/noLateral/runs/',filespec,'_i1/',filespec,'_i1_output',yearspec];
% configfile2 = ['../runs/noLateral/runs/',filespec,'_i1/',filespec,'_i1_settings'];

%% Non_Lake
load(outputfile1); load(configfile1);
txt_setting=['Infil: ',num2str(PARA.modules.infiltration),' Xice: ',num2str(PARA.modules.xice),' heat ex: ',num2str(PARA.modules.exchange_heat),' water ex: ',num2str(PARA.modules.exchange_water),' snow ex: ',num2str(PARA.modules.exchange_snow), ' LakeDepth: ',num2str(PARA.water.depth)];
txt_forcing=[PARA.forcing.filename(1:end-4),' (',num2str(PARA.forcing.rain_fraction),' ',num2str(PARA.forcing.snow_fraction),')'];

%PARA=p1; GRID=g1; OUT=o1; 
%z_soil=GRID.soil.soilGrid; ub_soil=GRID.soil.cT_domain_ub;
Soil_NL_ub = GRID.soil.cT_domain_ub; Soil_NL_lb = GRID.soil.cT_domain_lb;
zgeneral_NL=GRID.general.cT_grid;
zsoil_NL = GRID.general.cT_grid(Soil_NL_ub:Soil_NL_lb);
[~, z10m_NL]=min(abs(zgeneral_NL-10)); [~, z15m_NL]=min(abs(zgeneral_NL-15)); [~, z20m_NL]=min(abs(zgeneral_NL-20));
dl=10; % number of lines to plot
dz= round(length(zgeneral_NL)/dl);
tstamp=OUT.timestamp;
T_NL=OUT.cryoGrid3; T_NLtm = mean(T_NL,2);
TSoil_NL=T_NL(Soil_NL_ub:Soil_NL_lb,:); 
%TAir=T_NL(GRID.air.cT_domain_ub:GRID.air.cT_domain_lb,:); 
TSnow=T_NL(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb,:);
TSnow_lb=T_NL(GRID.lake.water.cT_domain_ub-1,:);
%OUT.soil.topPosition
% Forcing
tspan=FORCING.data.t_span; [ ~ ,idx1 ] = min( abs( tstamp(1) - tspan )); [ ~ ,idx2 ] = min( abs( tstamp(end) - tspan ));
tstamp_forcing=FORCING.data.t_span(idx1:idx2);
Tair=FORCING.data.Tair(idx1:idx2);

Qlat_NL=OUT.EB.Q_lateral; Qlat_NLtm = mean(Qlat_NL,2);

%%
co=cool(dl); % new Color Order
set(groot,'defaultAxesColorOrder',co)

    figure(1)
plot(tstamp,T_NL)
grid on; datetick; xlabel('timestamp'); ylabel('T_{NL} (°C)'); grid on;
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')

    figure(2)
%plot(tstamp,TAir,'b',
%tstamp,TSnow,'y',tstamp,TLakeWater,'b',tstamp,TLakeIce,'c',tstamp,TSoil,'g');plot(tstamp,TSnow,'k--',tstamp,TSoil(1:dz:end,:),'g')
plot(tstamp,TSoil_NL(1:dz:end,:),'g')
%hold on
% if(~isempty(GRID.lake.water.cT_domain_ub)); plot(tstamp,TLakeWater(1:dz:end,:),'b'); end
% if(~isempty(GRID.lake.ice.cT_domain_ub)); plot(tstamp,TLakeIce(1:dz:end,:),'c'); end
% hold off
grid on; datetick; xlabel('timestamp'); ylabel('T Soil (°C)'); grid on;

%% lateral heat fluxes
figure(4)
colormap(cool)
plot(tstamp,Qlat_NL(Soil_NL_ub:dz:end,:))
%plot(tstamp,Qlat_NL(Soil_NL_ub:Soil_NL_ub+2,:)')
grid on; datetick; ylabel('Qlat NL (W/m2)')

% figure(5)
% ti=100; z1=83; z2=300;
%     subplot(1,2,1)
% plot(Qlat_NL(z1:z2,ti),zgeneral(z1:z2))
% hold on; plot(Qlat_NLtm(z1:z2),zgeneral(z1:z2),'--'); hold off 
% set(gca,'Ydir','reverse'); xlabel('Qlat Non-LAKE (W/m2)'); ylabel('z (m)'); grid on
%     subplot(1,2,2)
% plot(T_L(z1:z2,ti)-T_NL(z1:z2,ti),zgeneral(z1:z2))
% hold on; plot(T_Ltm(z1:z2)-T_NLtm(z1:z2),zgeneral(z1:z2),'--'); hold off 
% set(gca,'Ydir','reverse'); xlabel('T_L-T_{NL} (°C)'); ylabel('z (m)'); grid on
% 
% T_NL and T_L profiles
% figure(6)
%     subplot(1,2,1)
% plot(T_L(z1:z2,ti),zgeneral(z1:z2))
% hold on; plot(T_Ltm(z1:z2),zgeneral(z1:z2),'--'); hold off 
% set(gca,'Ydir','reverse'); xlabel('T_L (°C)'); ylabel('z (m)'); grid on
% axis([-10 10 -inf inf])
%    subplot(1,2,2)
% plot(T_NL(z1:z2,ti),zgeneral(z1:z2))
% hold on; plot(T_NLtm(z1:z2),zgeneral(z1:z2),'--'); hold off 
% set(gca,'Ydir','reverse'); xlabel('T_{NL} (°C)'); ylabel('z (m)'); grid on
% axis([-10 10 -inf inf])


%%  LAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(outputfile2); load(configfile2);
txt_setting=['Infil: ',num2str(PARA.modules.infiltration),' Xice: ',num2str(PARA.modules.xice),' heat ex: ',num2str(PARA.modules.exchange_heat),' water ex: ',num2str(PARA.modules.exchange_water),' snow ex: ',num2str(PARA.modules.exchange_snow), ' LakeDepth: ',num2str(PARA.water.depth)];
txt_forcing=[PARA.forcing.filename(1:end-4),' (',num2str(PARA.forcing.rain_fraction),' ',num2str(PARA.forcing.snow_fraction),')'];

%PARA=p2; GRID=g2; OUT=o2; 
n_LakeLayers=nansum(GRID.lake.water.cT_domain)+nansum(GRID.lake.ice.cT_domain);
Lake_ub = GRID.soil.cT_domain_ub-n_LakeLayers;
Lake_lb = GRID.soil.cT_domain_ub-1;
Soil_L_ub = GRID.soil.cT_domain_ub;
zgeneral_L=GRID.general.cT_grid;
zsoil = GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
[~, z10m_L]=min(abs(zgeneral_L-10)); [~, z15m_L]=min(abs(zgeneral_L-15)); [~, z20m_L]=min(abs(zgeneral_L-20));
tstamp=OUT.timestamp;

T_L=OUT.cryoGrid3; T_Ltm = mean(T_L,2);
TSoil_L=T_L(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:); %TLakeWater=T_L(GRID.lake.water.cT_domain_ub:GRID.lake.water.cT_domain_lb,:);  TLakeIce=T_L(GRID.lake.ice.cT_domain_ub:GRID.lake.ice.cT_domain_lb,:); 
%TAir=T_L(GRID.air.cT_domain_ub:GRID.air.cT_domain_lb,:); 
TSnow=T_L(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb,:);
TSnow_lb=T_L(GRID.lake.water.cT_domain_ub-1,:);
TLake = T_L(Lake_ub:Lake_lb,:); 
TLake_ub = T_L(GRID.soil.cT_domain_ub-n_LakeLayers,:); % uppermost lake layer
TLake_lb = T_L(GRID.soil.cT_domain_ub-1,:); % lowermost lake layer
TaboveLake = T_L(GRID.soil.cT_domain_ub-(n_LakeLayers+1),:); 
TSoil_L_ub = T_L(GRID.soil.cT_domain_ub,:); 

Qlat_L=OUT.EB.Q_lateral; Qlat_Ltm = mean(Qlat_L,2);

%%
    figure(11)
plot(tstamp,T_L)
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); grid on;
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')

    figure(12) % T Soil and T Lake
%plot(tstamp,TAir,'b', tstamp,TSnow,'y',tstamp,TLakeWater,'b',tstamp,TLakeIce,'c',tstamp,TSoil_L,'g')
%plot(tstamp,TSnow,'k--',tstamp,TSoil_L(1:dz:end,:),'g')
plot(tstamp,TSoil_L(1:dz:end,:),'g')
hold on
if(~isempty(GRID.lake.water.cT_domain_ub)); plot(tstamp,TLake(1:dz:end,:),'b'); end
%if(~isempty(GRID.lake.water.cT_domain_ub)); plot(tstamp,TLakeWater(1:dz:end,:),'b'); end
%if(~isempty(GRID.lake.ice.cT_domain_ub)); plot(tstamp,TLakeIce(1:dz:end,:),'c'); end
hold off
grid on; datetick; xlabel('timestamp'); ylabel('T (soil & lake)(°C)'); 
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')
 
figure(13) % T Lake
%plot(tstamp,TLake(1:49:end,:),'b',tstamp,TaboveLake,'g',tstamp,TSoil_L_ub,'r')
plot(tstamp,TLake_ub,'b',tstamp,TLake_lb,'c',tstamp,TaboveLake,'mx',tstamp,TSoil_L_ub,'g')
legend('Tlake ub','Tlake lb','TaboveLake','TSoil ub')
hold on; plot(tstamp,TLake(1:10:end,:),'b--'); hold off
grid on; datetick; xlabel('timestamp'); ylabel('T (lake + ...)(°C)'); 
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')


%% lateral heat fluxes
ti=100; 

figure(14)
plot(tstamp,Qlat_L(Lake_ub:dz:end,:)')
grid on; datetick; ylabel('Qlat LAKE (W/m2)')
legend('1','2','3','4','5','6','7','8','9','10')

figure(15)
    subplot(1,2,1)
plot(Qlat_L(Lake_ub:end,ti),zgeneral_L(Lake_ub:end),'b--',Qlat_Ltm(Lake_ub:end),zgeneral_L(Lake_ub:end),'b');
legend('ti','tm')
set(gca,'Ydir','reverse'); xlabel('Qlat LAKE (W/m2)'); ylabel('z (m)'); grid on; title([filespec,' ',yearspec])
    subplot(1,2,2)
% plot(T_L(z1:z2,ti)-T_NL(z1:z2,ti),zgeneral(z1:z2),'b--',T_Ltm(z1:z2)-T_NLtm(z1:z2),zgeneral(z1:z2),'b'    is not on same depth!!
%legend('ti','tm')
%set(gca,'Ydir','reverse'); xlabel('T_L-T_{NL} (°C)'); ylabel('z (m)'); grid on

% T_NL and T_L profiles
figure(16)
    subplot(1,3,1)
plot(T_L(Lake_ub:end,ti),zgeneral_L(Lake_ub:end),'b--',T_Ltm(Lake_ub:end),zgeneral_L(Lake_ub:end),'b',T_NL(Soil_NL_ub:end,ti),zgeneral_NL(Soil_NL_ub:end),'g--',T_NLtm(Soil_NL_ub:end),zgeneral_NL(Soil_NL_ub:end),'g')
legend('L(ti)','L(tm)','NL(ti)','NL(tm)')
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on; title([filespec,' ',yearspec])
axis([-12 10 -inf inf])
   subplot(1,3,2)
plot(T_L(Lake_ub:z20m_L,ti),zgeneral_L(Lake_ub:z20m_L),'b--',T_Ltm(Lake_ub:z20m_L),zgeneral_L(Lake_ub:z20m_L),'b',T_NL(Soil_NL_ub:z20m_NL,ti),zgeneral_NL(Soil_NL_ub:z20m_NL),'g--',T_NLtm(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'g')
legend('L(ti)','L(tm)','NL(ti)','NL(tm)')
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
   subplot(1,3,3)
plot(Qlat_L(Lake_ub:z20m_L,ti),zgeneral_L(Lake_ub:z20m_L),'b--',Qlat_Ltm(Lake_ub:z20m_L),zgeneral_L(Lake_ub:z20m_L),'b',-Qlat_NL(Soil_NL_ub:z20m_NL,ti),zgeneral_NL(Soil_NL_ub:z20m_NL),'g--',-Qlat_NLtm(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'g')
legend('L(ti)','L(tm)','-NL(ti)','-NL(tm)')
set(gca,'Ydir','reverse'); xlabel('Qlat (W/m2)'); ylabel('z (m)'); grid on


figure(17) % uppermost soil layer vs. uppermost lake layer
%plot(tstamp,T_NL(GRID.soil.cT_domain_ub-10:GRID.soil.cT_domain_ub,:),'g',tstamp,TLake_ub,'b')
plot(tstamp,TSoil_NL(1,:),'g',tstamp,T_L(Lake_ub,:),'b')
legend('top T(NL)','top T(lake)'); grid on; ylabel('T (°C)'); datetick

figure(21) % check whether Qlat(L) = - Qlat(NL)
plot(Qlat_NL(Soil_NL_ub:Soil_NL_lb,ti)+Qlat_L(Lake_ub:Soil_NL_lb,ti),zgeneral_L(Lake_ub:Soil_NL_lb),'b--',Qlat_NLtm(Soil_NL_ub:Soil_NL_lb)+Qlat_Ltm(Lake_ub:Soil_NL_lb),zgeneral_L(Lake_ub:Soil_NL_lb),'b')
legend('ti','tm')
set(gca,'Ydir','reverse'); xlabel('Qlat LAKE + Qlat NL (W/m2)'); ylabel('z (m)'); grid on; %title([filespec,' ',yearspec])
title(['Diff not on same level!!!!  z integrated Qlat L+NL at ti ',num2str(nansum(Qlat_L(:,ti)) + nansum(Qlat_NL(:,ti)) ),' W/m2'])
