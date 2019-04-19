% T_NL and T_L profiles

load_data=1; % load data which were generated on cluster

if load_data == 1 % zzz ... update data on cluster....
    load Data_Tz
else
    Tini_ssw=[0 1 2 10 20 100 2000; 10 -3 -4 -7 -10 -10 10; 6 6 0 -5 -10 -10 10];
    Tini_msw=[0 5 8 20 100 2000; 10 -3 -6 -10 -10 10; 5 5 0 -9 -10 10];
    Tini=Tini_ssw';

    yearspec='2012';

    LD=[1 5];
    LR=[10 100]; 
    LF=[10 25 50];

    n=0; % run index loop
    for i=1:length(LD) % loop over exps
    for j=1:length(LR)
    for k=1:length(LF)
        n=n+1;
        expspec_NL=['WD0_LD',num2str(LD(i)),'_LR',num2str(LR(j)),'_LF',num2str(LF(k))]; % experiment specifier non Lake
        expspec_L=['WD',num2str(LD(i)),'_LD',num2str(LD(i)),'_LR',num2str(LR(j)),'_LF',num2str(LF(k))];
    % experiment specifier lake
        % NL
        outputfile_NL = ['../runs/',expspec_NL,'/',expspec_NL,'_',yearspec];
        outputfile_NL_nolat = ['../../CryoGrid3_noLateral/runs/',expspec_NL,'/',expspec_NL,'_',yearspec];
        configfile_NL = ['../runs/',expspec_NL,'/',expspec_NL,'_settings'];
        % Lake
        outputfile_L = ['../runs/',expspec_L,'/',expspec_L,'_',yearspec];
        outputfile_L_nolat = ['../../CryoGrid3_noLateral/runs/',expspec_L,'/',expspec_L,'_',yearspec];
        configfile_L = ['../runs/',expspec_L,'/',expspec_L,'_settings'];
    % 1.) z levels
        if(n==1)
            load(configfile_NL); % NL
            Soil_NL_ub = GRID.soil.cT_domain_ub; Soil_NL_lb = GRID.soil.cT_domain_lb;
            zgeneral_NL=GRID.general.cT_grid; zsoil_NL = GRID.general.cT_grid(Soil_NL_ub:Soil_NL_lb);
            [~, z20m_NL]=min(abs(zgeneral_NL-30));

            load(configfile_L); % Lake 1m
            n_LakeLayers=nansum(GRID.lake.water.cT_domain)+nansum(GRID.lake.ice.cT_domain);
            Lake_1m_ub = GRID.soil.cT_domain_ub-n_LakeLayers; Lake_1m_lb = GRID.soil.cT_domain_ub-1;
            zgeneral_L1m=GRID.general.cT_grid;
            [~, z20m_L1m]=min(abs(zgeneral_L1m-30));
        elseif(n==7)
            load(configfile_L); % Lake 5m
            n_LakeLayers=nansum(GRID.lake.water.cT_domain)+nansum(GRID.lake.ice.cT_domain);
            Lake_5m_ub = GRID.soil.cT_domain_ub-n_LakeLayers; Lake_5m_lb = GRID.soil.cT_domain_ub-1;
            zgeneral_L5m=GRID.general.cT_grid;
            [~, z20m_L5m]=min(abs(zgeneral_L5m-30));
        end
    % 2.) generate T data
    % non Lake
        load(outputfile_NL); 
        T_NL=OUT.cryoGrid3; T_NLtm = mean(T_NL,2);
        eval(['T_',expspec_NL,'=T_NLtm;']);
       % no lateral heat flux    
        load(outputfile_NL_nolat);  
        T_NL=OUT.cryoGrid3; T_NLtm = mean(T_NL,2);
        eval(['T_',expspec_NL,'_nolat=T_NLtm;']);
    % Lake
        load(outputfile_L);
        T_L=OUT.cryoGrid3; T_Ltm = mean(T_L,2);
        eval(['T_',expspec_L,'=T_Ltm;']);
       % no lateral heat flux
        load(outputfile_L_nolat);
        T_L=OUT.cryoGrid3; T_Ltm = mean(T_L,2);
        eval(['T_',expspec_L,'_nolat=T_Ltm;']);
    end
    end
end 
save Data_Tz T* z* Lake_* yearspec expspec* Soil_NL_ub L*
end 

%%
figure(1) % LD 1m
    subplot(1,2,1)
%  SSW (LR 10m)
hold on
% LF 0.25
plot(T_WD0_LD1_LR10_LF25(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR10_LF25(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b')
% LF 0.1
plot(T_WD0_LD1_LR10_LF10(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k:',T_WD1_LD1_LR10_LF10(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b:') % WD0 is non Lake!
% LF 0.5
plot(T_WD0_LD1_LR10_LF50(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k--',T_WD1_LD1_LR10_LF50(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b--')
% no lateral exchange
%plot(T_WD0_LD1_LR10_LF25_nolat(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR10_LF25_nolat(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'c')
% initial T 
%plot(Tini(1:5,2),Tini(1:5,1),'g*',Tini(1:5,3),Tini(1:5,1),'b*')
plot([-20 10],[LD(1) LD(1)],'b:','LineWidth',0.8) % lake depth level
hold off
legend('Soil','Lake','Location','SouthWest');  title(['lake radius ',num2str(LR(1)),'m'])  % title(['LD',num2str(LD(1)),' LR',num2str(LR(1))])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-20 10 0 10])

    subplot(1,2,2)
%  MSW (LR 100m)
hold on
plot(T_WD0_LD1_LR100_LF25(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR100_LF25(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b')
plot(T_WD0_LD1_LR100_LF10(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k:',T_WD1_LD1_LR100_LF10(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b:')
plot(T_WD0_LD1_LR100_LF50(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k--',T_WD1_LD1_LR100_LF50(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b--')
% no lateral exchange
%plot(T_WD0_LD1_LR100_LF25_nolat(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR100_LF25_nolat(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'c')
%plot(Tini(1:5,2),Tini(1:5,1),'g*',Tini(1:5,3),Tini(1:5,1),'b*')
plot([-20 10],[LD(1) LD(1)],'b:','LineWidth',0.8) % lake depth level
hold off
legend('Soil','Lake','Location','SouthWest');  title(['lake radius ',num2str(LR(2)),'m']) % title(['LD',num2str(LD(1)),' LR',num2str(LR(2))])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-20 10 0 10])

print -dtiff Tz_LD1m


%%
figure(2) % LD 5m
    subplot(1,2,1)
%  SSW (LR 10m)
hold on
% LF 0.25
plot(T_WD0_LD5_LR10_LF25(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD5_LD5_LR10_LF25(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'b')
% LF 0.1
plot(T_WD0_LD5_LR10_LF10(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k:',T_WD5_LD5_LR10_LF10(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'b:') % WD0 is non Lake!
% LF 0.5
plot(T_WD0_LD5_LR10_LF50(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k--',T_WD5_LD5_LR10_LF50(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'b--')
% no lateral exchange
%plot(T_WD0_LD5_LR10_LF25_nolat(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD5_LD5_LR10_LF25_nolat(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'c')
% initial T 
%plot(Tini(1:5,2),Tini(1:5,1),'g*',Tini(1:5,3),Tini(1:5,1),'b*')
plot([-20 10],[LD(2) LD(2)],'b:','LineWidth',0.8) % lake depth level
hold off
legend('Soil','Lake','Location','SouthWest');  title(['lake radius ',num2str(LR(1)),'m'])  % title(['LD',num2str(LD(1)),' LR',num2str(LR(1))])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-20 10 0 10])

    subplot(1,2,2)
%  MSW (LR 100m)
hold on
plot(T_WD0_LD5_LR100_LF25(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD5_LD5_LR100_LF25(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'b')
plot(T_WD0_LD5_LR100_LF10(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k:',T_WD5_LD5_LR100_LF10(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'b:')
plot(T_WD0_LD5_LR100_LF50(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k--',T_WD5_LD5_LR100_LF50(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'b--')
% no lateral exchange
%plot(T_WD0_LD5_LR100_LF25_nolat(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR100_LF25_nolat(Lake_5m_ub:z20m_L5m),zgeneral_L5m(Lake_5m_ub:z20m_L5m),'c')
%plot(Tini(1:5,2),Tini(1:5,1),'g*',Tini(1:5,3),Tini(1:5,1),'b*')
plot([-20 10],[LD(2) LD(2)],'b:','LineWidth',0.8) % lake depth level
hold off
legend('Soil','Lake','Location','SouthWest');  title(['lake radius ',num2str(LR(2)),'m']) % title(['LD',num2str(LD(1)),' LR',num2str(LR(2))])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-20 10 0 10])

print -dtiff Tz_LD5m


%%
figure(3) % LD 1m
    subplot(1,2,1)
%  SSW (LR 10m)
hold on
% LF 0.25
plot(T_WD0_LD1_LR10_LF25(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR10_LF25(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b')
% no lateral exchange
plot(T_WD0_LD1_LR10_LF25_nolat(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k--',T_WD1_LD1_LR10_LF25_nolat(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b--')
plot([-20 10],[LD(1) LD(1)],'b:','LineWidth',0.8) % lake depth level
hold off
legend('Soil','Lake','Location','SouthWest');  title(['lake radius ',num2str(LR(1)),'m'])  % title(['LD',num2str(LD(1)),' LR',num2str(LR(1))])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-20 10 0 10])

    subplot(1,2,2)
%  MSW (LR 100m)
hold on
plot(T_WD0_LD1_LR100_LF25(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k',T_WD1_LD1_LR100_LF25(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b')
% no lateral exchange
plot(T_WD0_LD1_LR100_LF25_nolat(Soil_NL_ub:z20m_NL),zgeneral_NL(Soil_NL_ub:z20m_NL),'k--',T_WD1_LD1_LR100_LF25_nolat(Lake_1m_ub:z20m_L1m),zgeneral_L1m(Lake_1m_ub:z20m_L1m),'b--')
plot([-20 10],[LD(1) LD(1)],'b:','LineWidth',0.8) % lake depth level
hold off
legend('Soil','Lake','Location','SouthWest');  title(['lake radius ',num2str(LR(2)),'m']) % title(['LD',num2str(LD(1)),' LR',num2str(LR(2))])
set(gca,'Ydir','reverse'); xlabel('T (°C)'); ylabel('z (m)'); grid on
axis([-20 10 0 10])

print -dtiff Tz_LD1m_nolat
