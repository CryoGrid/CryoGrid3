% T_NL and T_L profiles

%filespec='MSW_xH1_xW0_xS0_infil1_xice0_rF1_sF1';
%filespec_b='MSW_xH0_xW0_xS0_infil1_xice0_rF1_sF1'; % no lateral heat exchange
filespec='xH1';
filespec_b='xH0';
yearspec='1971';

Tini_ssw=[0 1 2 10 20 100 2000; 10 -3 -4 -7 -10 -10 10; 6 6 0 -5 -10 -10 10];
Tini_msw=[0 5 8 20 100 2000; 10 -3 -6 -10 -10 10; 5 5 0 -9 -10 10];
Tini=Tini_ssw';

%% Non_Lake
% with lateral exchange
%outputfile1 = ['../runs/',filespec,'_i1/',filespec,'_i1_realization1_output',yearspec]; configfile1 = ['../runs/',filespec,'_i1/',filespec,'_i1_settings'];
outputfile1 = ['../runs/WD0_',filespec,'/WD0_',filespec,'_output',yearspec];  configfile1 = ['../runs/WD0_',filespec,'/WD0_',filespec,'_settings'];
load(outputfile1); load(configfile1);
Soil_NL_ub = GRID.soil.cT_domain_ub; Soil_NL_lb = GRID.soil.cT_domain_lb;
zgeneral_NL=GRID.general.cT_grid; zsoil_NL = GRID.general.cT_grid(Soil_NL_ub:Soil_NL_lb);
[~, z10m_NL]=min(abs(zgeneral_NL-10)); [~, z15m_NL]=min(abs(zgeneral_NL-15)); [~, z20m_NL]=min(abs(zgeneral_NL-20));

T_NL=OUT.cryoGrid3; T_NLtm = mean(T_NL,2);

% wo lateral exchange
%outputfile1b = ['../runs/noLateral/runs/',filespec_b,'_i1/',filespec_b,'_i1_realization1_output',yearspec]; configfile1b = ['../runs/noLateral/runs/',filespec_b,'_i1/',filespec_b,'_i1_settings'];
outputfile1b = ['../runs/noLateral/runs/WD0_',filespec_b,'/WD0_',filespec_b,'_output',yearspec]; configfile1b = ['../runs/noLateral/runs/WD0_',filespec_b,'/WD0_',filespec_b,'_settings'];
load(outputfile1b); load(configfile1b);

T_NL_xH0=OUT.cryoGrid3; T_NLtm_xH0 = mean(T_NL_xH0,2); % no lateral heat exchange

%%  LAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputfile2 = ['../runs/',filespec,'_i2/',filespec,'_i2_realization2_output',yearspec]; configfile2 = ['../runs/',filespec,'_i2/',filespec,'_i2_settings'];
outputfile2 = ['../runs/WD1_',filespec,'/WD1_',filespec,'_output',yearspec];  configfile2 = ['../runs/WD1_',filespec,'/WD1_',filespec,'_settings'];
load(outputfile2); load(configfile2)

n_LakeLayers=nansum(GRID.lake.water.cT_domain)+nansum(GRID.lake.ice.cT_domain);
Lake_ub = GRID.soil.cT_domain_ub-n_LakeLayers; Lake_lb = GRID.soil.cT_domain_ub-1;
Soil_L_ub = GRID.soil.cT_domain_ub;
zgeneral_L=GRID.general.cT_grid;
zsoil = GRID.general.cT_grid(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_lb);
[~, z10m_L]=min(abs(zgeneral_L-10)); [~, z15m_L]=min(abs(zgeneral_L-15)); [~, z20m_L]=min(abs(zgeneral_L-20));

T_L=OUT.cryoGrid3; T_Ltm = mean(T_L,2);

% wo lateral exchange
outputfile2b = ['../runs/noLateral/runs/WD1_',filespec_b,'/WD1_',filespec_b,'_output',yearspec]; configfile2b = ['../runs/noLateral/runs/WD1_',filespec_b,'/WD1_',filespec_b,'_settings'];
load(outputfile2b); load(configfile2b);

T_L_xH0=OUT.cryoGrid3; T_Ltm_xH0 = mean(T_L_xH0,2);

%%
ti=100;
figure
    subplot(1,2,1)
plot(T_L(Lake_ub:z20m_L,ti),zgeneral_L(Lake_ub:z20m_L),'b',T_NL(Soil_NL_ub:z20m_NL,ti),zgeneral_NL(Soil_NL_ub:z20m_NL),'g', ...
     T_L_xH0(Lake_ub:z20m_L,ti),zgeneral_L(Lake_ub:z20m_L),'b--',T_NL_xH0(Soil_NL_ub:z20m_NL,ti),zgeneral_NL(Soil_NL_ub:z20m_NL),'g--') 
hold on; plot(Tini(1:5,2),Tini(1:5,1),'g*',Tini(1:5,3),Tini(1:5,1),'b*'); hold off
 legend('L','NL','Location','SouthWest'); title(['T(z) at ti  ',num2str(filespec(1:3)),'  ',num2str(yearspec)])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-30 10 -inf inf])
    subplot(1,2,2)
plot(T_Ltm(Lake_ub:z20m_L),zgeneral_L(Lake_ub:z20m_L),'b',T_NLtm(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'g', ...
     T_Ltm_xH0(Lake_ub:z20m_L),zgeneral_L(Lake_ub:z20m_L),'b--',T_NLtm_xH0(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'g--')
hold on; plot(Tini(1:5,2),Tini(1:5,1),'g*',Tini(1:5,3),Tini(1:5,1),'b*'); hold off
legend('L','NL','Location','SouthWest');  title(['T(z) annual mean  ',num2str(filespec(1:3)),'  ',num2str(yearspec)])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-30 10 -inf inf])
