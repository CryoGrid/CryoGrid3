% plot T timeseries for various depth

% OUT=OUTS{2}.OUT; GRID=OUTS{2}.GRID; PARA=OUTS{2}.PARA

%p1=PARA{1}; p2=PARA{2}; g1=GRID{1}; g2=GRID{2}; o1=OUT{1}; o2=OUT{2};

dz=10; % increment used to plot vertical levels
%% Non_Lake

% no infil no ex no lake
% outputfile1 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i1/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i1_realization1_output1980.mat'
% configfile1 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i1/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i1_realization1_settings.mat'
% outputfile2 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i2/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i2_realization2_output1980.mat'
% configfile2 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i2/LAKE-MPI_xH1_xW0_xS0_infil0_xice0_rF1_sF1_i2_realization2_settings.mat'
% 
outputfile1 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i1/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i1_realization1_output1980.mat'
configfile1 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i1/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i1_realization1_settings.mat'
outputfile2 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i2/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i2_realization2_output1980.mat'
configfile2 = '../runs/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i2/LAKE-MPI_xH1_xW0_xS0_infil1_xice0_rF1_sF1_i2_realization2_settings.mat'

load(outputfile1); load(configfile1);

%PARA=p1; GRID=g1; OUT=o1; 
%z_soil=GRID.soil.soilGrid; ub_soil=GRID.soil.cT_domain_ub;
zsoil = GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
tstamp=OUT.timestamp;
T=OUT.cryoGrid3; TSoil=T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:); %TLakeWater=T(GRID.lake.water.cT_domain_ub:GRID.lake.water.cT_domain_lb,:);  TLakeIce=T(GRID.lake.ice.cT_domain_ub:GRID.lake.ice.cT_domain_lb,:); 
%OUT.soil.topPosition

%TAir=T(GRID.air.cT_domain_ub:GRID.air.cT_domain_lb,:); 

TSnow=T(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb,:);
TSnow_lb=T(GRID.lake.water.cT_domain_ub-1,:);

% Forcing
 tspan=FORCING.data.t_span;
[ ~ ,idx1 ] = min( abs( tstamp(1) - tspan )); [ ~ ,idx2 ] = min( abs( tstamp(end) - tspan ));
tstamp_forcing=FORCING.data.t_span(idx1:idx2);
Tair=FORCING.data.Tair(idx1:idx2);

txt_setting=['Infil: ',num2str(PARA.modules.infiltration),' Xice: ',num2str(PARA.modules.xice),' heat ex: ',num2str(PARA.modules.exchange_heat),' water ex: ',num2str(PARA.modules.exchange_water),' snow ex: ',num2str(PARA.modules.exchange_snow), ' LakeDepth: ',num2str(PARA.water.depth)];
txt_forcing=[PARA.forcing.filename(1:end-4),' (',num2str(PARA.forcing.rain_fraction),' ',num2str(PARA.forcing.snow_fraction),')'];

    figure(1)
plot(tstamp,T)
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); grid on;
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')

    figure(2)
%plot(tstamp,TAir,'b', tstamp,TSnow,'y',tstamp,TLakeWater,'b',tstamp,TLakeIce,'c',tstamp,TSoil,'g')
%plot(tstamp,TSnow,'k--',tstamp,TSoil(1:dz:end,:),'g')
plot(tstamp,TSoil(1:dz:end,:),'g')
%hold on
% if(~isempty(GRID.lake.water.cT_domain_ub)); plot(tstamp,TLakeWater(1:dz:end,:),'b'); end
% if(~isempty(GRID.lake.ice.cT_domain_ub)); plot(tstamp,TLakeIce(1:dz:end,:),'c'); end
% hold off
grid on; datetick; xlabel('timestamp'); ylabel('T Soil (°C)'); grid on;

%%  LAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARA=p2; GRID=g2; OUT=o2; 
load(outputfile2); load(configfile2);
zsoil = GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
tstamp=OUT.timestamp;
T=OUT.cryoGrid3; TSoil=T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:); %TLakeWater=T(GRID.lake.water.cT_domain_ub:GRID.lake.water.cT_domain_lb,:);  TLakeIce=T(GRID.lake.ice.cT_domain_ub:GRID.lake.ice.cT_domain_lb,:); 
%TAir=T(GRID.air.cT_domain_ub:GRID.air.cT_domain_lb,:); 
TSnow=T(GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb,:);
TSnow_lb=T(GRID.lake.water.cT_domain_ub-1,:);
TLake = T(GRID.soil.cT_domain_ub-50:GRID.soil.cT_domain_ub-1,:); % lake 50 layers
Tx1 = T(GRID.soil.cT_domain_ub-51,:); 
Tx2 = T(GRID.soil.cT_domain_ub+1,:); 

    figure(11)
txt_setting=['Infil: ',num2str(PARA.modules.infiltration),' Xice: ',num2str(PARA.modules.xice),' heat ex: ',num2str(PARA.modules.exchange_heat),' water ex: ',num2str(PARA.modules.exchange_water),' snow ex: ',num2str(PARA.modules.exchange_snow), ' LakeDepth: ',num2str(PARA.water.depth)];
txt_forcing=[PARA.forcing.filename(1:end-4),' (',num2str(PARA.forcing.rain_fraction),' ',num2str(PARA.forcing.snow_fraction),')'];
plot(tstamp,T)
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); grid on;
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')

    figure(12)
%plot(tstamp,TAir,'b', tstamp,TSnow,'y',tstamp,TLakeWater,'b',tstamp,TLakeIce,'c',tstamp,TSoil,'g')
%plot(tstamp,TSnow,'k--',tstamp,TSoil(1:dz:end,:),'g')
plot(tstamp,TSoil(1:dz:end,:),'g')
hold on
if(~isempty(GRID.lake.water.cT_domain_ub)); plot(tstamp,TLake(1:dz:end,:),'b'); end
%if(~isempty(GRID.lake.water.cT_domain_ub)); plot(tstamp,TLakeWater(1:dz:end,:),'b'); end
%if(~isempty(GRID.lake.ice.cT_domain_ub)); plot(tstamp,TLakeIce(1:dz:end,:),'c'); end
hold off
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); 
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')
 
%%
%[ ~ ,idx ] = min( abs( datenum(1979, 7, 1 ) - OUT.timestamp() )) % find index of date
    figure(30)
year=1980;
[~,tJan]=min(abs(datenum(year,1,1)-tstamp())); [~,tFeb]=min(abs(datenum(year,2,1)-tstamp())); [~,tMar]=min(abs(datenum(year,3,1)-tstamp())); [~,tApr]=min(abs(datenum(year,4,1)-tstamp())); [~,tMay]=min(abs(datenum(year,5,1)-tstamp())); [~,tJun]=min(abs(datenum(year,6,1)-tstamp())); 
[~,tJul]=min(abs(datenum(year,7,1)-tstamp())); [~,tAug]=min(abs(datenum(year,8,1)-tstamp())); [~,tSep]=min(abs(datenum(year,9,1)-tstamp())); [~,tOct]=min(abs(datenum(year,4,10)-tstamp())); [~,tNov]=min(abs(datenum(year,5,11)-tstamp())); [~,tDec]=min(abs(datenum(year,12,1)-tstamp())); 
 
z2m=101; z10m=261; % corresponds to 2m, 10m
zend=z2m;
plot(TSoil(1:zend,tJan),zsoil(1:zend),TSoil(1:zend,tFeb),zsoil(1:zend),TSoil(1:zend,tMar),zsoil(1:zend),TSoil(1:zend,tApr),zsoil(1:zend),TSoil(1:zend,tMay),zsoil(1:zend),TSoil(1:zend,tJun),zsoil(1:zend),...
     TSoil(1:zend,tJul),zsoil(1:zend),TSoil(1:zend,tAug),zsoil(1:zend),TSoil(1:zend,tSep),zsoil(1:zend),TSoil(1:zend,tOct),zsoil(1:zend),TSoil(1:zend,tNov),zsoil(1:zend),TSoil(1:zend,tDec),zsoil(1:zend))
hold on; 
set(gca,'Ydir','reverse'); xlabel('T'); ylabel('z'); grid on
legend(datestr(tstamp(tJan)),datestr(tstamp(tFeb)),datestr(tstamp(tMar)),datestr(tstamp(tApr)),datestr(tstamp(tMay)),datestr(tstamp(tJun)),datestr(tstamp(tJul)),datestr(tstamp(tAug)),datestr(tstamp(tSep)),datestr(tstamp(tOct)),datestr(tstamp(tNov)),datestr(tstamp(tDec)))


figure(99)
hold on
plot(tstamp,TLake(1:49:end,:),'b',tstamp,Tx1,'g',tstamp,Tx2,'r')
grid on; datetick
hold off
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); 
xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
title(txt_forcing, 'Interpreter', 'none')