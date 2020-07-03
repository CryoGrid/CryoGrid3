% plot time series of frost and tables
% first run extract_timeseries*.m or load Data and execute plot_FT_WT.m to get zsoil and idx10m, idx20
% assume that no change in soil levels!

numTiles=5;
%dataDir='Runs_paper'
%dataDir='Data5'
dataDir='Data/Runs_paper_final'
%dataDir='Data/Runs_paper_final_DP'

SetLoc='Prudhoe'; 
%SetLoc='Drew'; %
ExpSet='Drain0cmPond5cm'
%ExpSet='Drain0cmPond5cm_snow30peat15_ERA'
%Scen='RCP85';
Scen='CTR'
%Scen='CTRERA';
loadData=1;    
plotSwitch=1;

if(loadData==1)
    load([dataDir,'/Data',num2str(numTiles),'_Tsoil_10m20m_',SetLoc,'_',Scen,'_',ExpSet])
    load([dataDir,'/Data',num2str(numTiles),'_FTWT_',SetLoc,'_',Scen,'_',ExpSet])
    load([dataDir,'/Data',num2str(numTiles),'_Q_',SetLoc,'_',Scen,'_',ExpSet])
    load([dataDir,'/Data',num2str(numTiles),'_F_',SetLoc,'_',Scen])
    %load([dataDir,'/Data5_Tsoil_10m20m_',SetLoc,'_',Scen,'_',ExpSet]); load([dataDir,'/Data5_FTWT_',SetLoc,'_',Scen,'_',ExpSet]); load([dataDir,'/Data5_F_',SetLoc,'_',Scen])
end
time = datetime(ts,'ConvertFrom','datenum');
tt1=datevec(ts(1)); tt2=datevec(ts(end)); 
%t1=1980; t2=1981;
t1=tt1(1); t2=tt2(1); % start and endyear
%t1=1980; t2=2050;

%% 
    figure(3) % ALD
hold on
%plot(t1:t2,FTmin_T1,'+-',t1:t2,FTmin_T2,'+-',t1:t2,FTmin_T3,'+-',t1:t2,FTmin_T4,'+-',t1:t2,FTmin_T5,'+-')
if(plotSwitch==1)
    if(numTiles==5)
        h1=plot(t1:t2,FTmin_T1,'c.',t1:t2,FTmin_T2,'b.',t1:t2,FTmin_T3,'r.',t1:t2,FTmin_T4,'g.',t1:t2,FTmin_T5,'k.','MarkerSize',20);
    elseif(numTiles==30)
        h1=plot(t1:t2,FTmin_T3,'cv',t1:t2,FTmin_T8,'bv',t1:t2,FTmin_T13,'rv',t1:t2,FTmin_T18,'gv',t1:t2,FTmin_T23,'kv','MarkerSize',2);        
    end
   % set(h1,'LineWidth',8)
elseif(plotSwitch==2)
    if(numTiles==5)
        plot(t1:t2,FTmin_T1,'c--',t1:t2,FTmin_T2,'b--',t1:t2,FTmin_T3,'r--',t1:t2,FTmin_T4,'g--',t1:t2,FTmin_T5,'k--')
    elseif(numTiles==30)
        plot(t1:t2,FTmin_T3,'c--',t1:t2,FTmin_T8,'b--',t1:t2,FTmin_T13,'r--',t1:t2,FTmin_T18,'g--',t1:t2,FTmin_T23,'k--')        
    end
end
%set(h1,'MarkerSize',10)
h2=plot([t1 t2],[-1.5 -1.5],'--k'); % embankment depth
set(h2,'LineWidth',1.5)
xlim([2000 2075]); ylim([-10 1]);
legend('T1','T2','T3','T4','T5','EB base','Location','SouthWest'); grid on; grid minor; ylabel('frost table [m]'); title([SetLoc,' ',Scen])
xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet,'Interpreter', 'none')
%plot([0 10],[-PARA.IS.EHbg -PARA.IS.EHbg],'--k','LineWidth',1.5) % embankment depth
hold off

    figure(33) % plot FTmin 5tiles vs 30tiles
hold on
if(numTiles==5)
    plot(t1:t2,FTmin_T1,'c.',t1:t2,FTmin_T2,'b.',t1:t2,FTmin_T3,'r.',t1:t2,FTmin_T4,'g.',t1:t2,FTmin_T5,'k.','MarkerSize',20);
elseif(numTiles==30)
    plot(t1:t2,FTmin_T3,'cv',t1:t2,FTmin_T8,'bv',t1:t2,FTmin_T13,'rv',t1:t2,FTmin_T18,'gv',t1:t2,FTmin_T23,'kv','MarkerSize',2);        
end
hold off

    figure(1) % Tsoil 10m and 20m
plot(time,Tsoil_10m_tiles,time,Tsoil_20m_tiles,'--',timeF,movmean(TairF,365*24),'m')
%xlim([datetime(1980,1,1) datetime(2020,1,1)]); ylim([-13 2])
legend('Tsoil 10m','Tsoil 20m','Tair')
ylabel('Tsoil_{10/20} and mean Tair');grid on; grid minor; title([SetLoc,' ',Scen])
xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet,'Interpreter', 'none')

    figure(2) % frost table
plot(time,FT_T1,time,FT_T2,time,FT_T3,time,FT_T4,time,FT_T5)
legend('T1','T2','T3','T4','T5'); grid on; grid minor; ylabel('frost table [m]'); title([SetLoc,' ',Scen])
xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet,'Interpreter', 'none')

    figure(4) % water table
plot(time,WT_T1,time,WT_T2,time,WT_T3,time,WT_T4,time,WT_T5)
legend('T1','T2','T3','T4','T5'); grid on; grid minor; ylabel('water table [m]'); title([SetLoc,' ',Scen])
xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet,'Interpreter', 'none')

    figure(5) % soil level
plot(time,SL_T1,time,SL_T2,time,SL_T3,time,SL_T4,time,SL_T5)
legend('T1','T2','T3','T4','T5'); grid on; grid minor; ylabel('soil level [m]'); title([SetLoc,' ',Scen])
xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet,'Interpreter', 'none')

    figure(6) % Tair forcing
plot(timeF,TairF,timeF,movmean(TairF,365*24))
%xlim([datetime(1980,1,1) datetime(2020,1,1)])
xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet,'Interpreter', 'none')
title([SetLoc,' ',Scen]); grid on; grid minor; ylabel('Tair forcing')

% %% 24 tiles
% if(loadData==1)
%     load(['Data24/Data24_Tsoil_10m20m_',SetLoc,'_',Scen,'_',ExpSet])
%     load(['Data24/Data24_FTWT_',SetLoc,'_',Scen,'_',ExpSet])
%     load(['Data24/Data24_F_',SetLoc,'_',Scen])
% end
% time = datetime(ts,'ConvertFrom','datenum');
% tt1=datevec(ts(1)); tt2=datevec(ts(end)); 
% t1=tt1(1); t2=tt2(1); % start and endyear
% 
% %%
%     figure(11) % Tsoil 10m and 20m
% hold on
% plot(time,Tsoil_10m_tiles,time,Tsoil_20m_tiles,'--',time_F,movmean(Tair_F,365*24),'m')
% xlim([datetime(1980,1,1) datetime(2017,1,1)])
% ylabel('Tsoil_{10/20} and mean Tair');grid on; grid minor; title([SetLoc,' ',Scen])
% xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet)
% 
%     figure(13) % ALD
% hold on
% %plot(t1:t2,FTmin_T1,'+-',t1:t2,FTmin_T2,'+-',t1:t2,FTmin_T3,'+-',t1:t2,FTmin_T4,'+-',t1:t2,FTmin_T5,'+-')
% plot(t1:t2,FTmin_T1,t1:t2,FTmin_T2,t1:t2,FTmin_T3,t1:t2,FTmin_T4,t1:t2,FTmin_T5,t1:t2,FTmin_T6,t1:t2,FTmin_T7,t1:t2,FTmin_T8,t1:t2,FTmin_T9,t1:t2,FTmin_T10,t1:t2,FTmin_T11,t1:t2,FTmin_T12)
% plot(t1:t2,FTmin_T13,t1:t2,FTmin_T14,t1:t2,FTmin_T15,t1:t2,FTmin_T16,t1:t2,FTmin_T17,t1:t2,FTmin_T18,t1:t2,FTmin_T19,t1:t2,FTmin_T20,t1:t2,FTmin_T21,t1:t2,FTmin_T22,t1:t2,FTmin_T23,t1:t2,FTmin_T24)
% legend('T1','T2','T3','T4','T5','EB base'); grid on; grid minor; ylabel('frost table [m]'); title([SetLoc,' ',Scen])
% xxx=xlim; yyy=ylim; text(xxx(1),0.8*yyy(2),ExpSet)
% 
