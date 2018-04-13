% plot T timeseries for various depth

% OUT=OUTS{2}.OUT; GRID=OUTS{2}.GRID; PARA=OUTS{2}.PARA

% GRID=g1; PARA=p1; OUT=o1;
%% 
figure(1)

outputfile = '../runs/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i1/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i1_realization1_output1980.mat'
configfile = '../runs/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i1/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i1_realization1_settings.mat'

load(outputfile); load(configfile);

%z_soil=GRID.soil.soilGrid; ub_soil=GRID.soil.cT_domain_ub;
zsoil = GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);

tstamp=OUT.timestamp;
T=OUT.cryoGrid3;
Tsoil=T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb,:);

txt_setting=['Infil: ',num2str(PARA.modules.infiltration),' Xice: ',num2str(PARA.modules.xice),' heat ex: ',num2str(PARA.modules.exchange_heat),' water ex: ',num2str(PARA.modules.exchange_water),' snow ex: ',num2str(PARA.modules.exchange_snow), ' LakeDepth: ',num2str(PARA.water.depth)];
txt_forcing=[PARA.forcing.filename(1:end-4),' (',num2str(PARA.forcing.rain_fraction),' ',num2str(PARA.forcing.snow_fraction),')'];
%plot(tstamp,T(ub_soil:end,:))
plot(tstamp,T(ub_soil:10:end,:))
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); grid on;

xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
%text(xxx(1),0.8*yyy(2),txt_forcing, 'Interpreter', 'none')
title(txt_forcing, 'Interpreter', 'none')

%%
figure(2)

outputfile = '../runs/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i2/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i2_realization2_output1980.mat'
configfile = '../runs/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i2/LAKE-MPI_xH0_xW0_xS0_infil0_xice0_rF1_sF1_i2_realization2_settings.mat'

load(outputfile); load(configfile);

%z_soil=GRID.soil.soilGrid; ub_soil=GRID.soil.cT_domain_ub;

tstamp=OUT.timestamp;
T=OUT.cryoGrid3;

txt_setting=['Infil: ',num2str(PARA.modules.infiltration),' Xice: ',num2str(PARA.modules.xice),' heat ex: ',num2str(PARA.modules.exchange_heat),' water ex: ',num2str(PARA.modules.exchange_water),' snow ex: ',num2str(PARA.modules.exchange_snow), ' LakeDepth: ',num2str(PARA.water.depth)];
txt_forcing=[PARA.forcing.filename(1:end-4),' (',num2str(PARA.forcing.rain_fraction),' ',num2str(PARA.forcing.snow_fraction),')'];
%plot(tstamp,T(ub_soil:end,:))
plot(tstamp,T(ub_soil:10:end,:))
grid on; datetick; xlabel('timestamp'); ylabel('T (°C)'); grid on;

xxx=xlim; yyy=ylim; text(xxx(1),0.9*yyy(2),txt_setting)
%text(xxx(1),0.8*yyy(2),txt_forcing, 'Interpreter', 'none')
title(txt_forcing, 'Interpreter', 'none')

%%
%[ ~ ,idx ] = min( abs( datenum(1979, 7, 1 ) - OUT.timestamp() )) % find index of date
figure(3)
year=1979;
[~,tJan]=min(abs(datenum(year,1,1)-OUT.timestamp())); [~,tApr]=min(abs(datenum(year,4,1)-OUT.timestamp()));  [~,tJul]=min(abs(datenum(year,7,1)-OUT.timestamp())); [~,tOct]=min(abs(datenum(year,10,1)-OUT.timestamp()));  %1. Jan 1. Apr  1.Jul  1.Oct
 
z2=101; z10=261; % corresponds to 2m, 10m
plot(Tsoil(1:z1,ti),zsoil(1:z1))
set(gca,'Ydir','reverse'); xlabel('T'); ylabel('z'); grid on
title(datestr(tstamp(ti)))
