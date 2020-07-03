% gravel road

OutDir = '.\Data\';
%%%OutTile = '2';
%Exp = ['GravelRoad_PrudoeBay_realization',OutTile,'_output'];
Exp = 'GravelRoad_PrudoeBay';
%% year x
%load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization2_output2000
%load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization2_finalState2000
OutYear = '2000'; OutTile='2';
OutputFile = [OutDir,Exp,'_realization',OutTile,'_output',OutYear]; FinalFile = [OutDir,Exp,'_realization',OutTile,'_finalState',OutYear];

%%%load Data\GravelRoad_PrudoeBay_realization2_output2000; load Data\GravelRoad_PrudoeBay_realization2_finalState2000
load(OutputFile); load(FinalFile)
ttt=datevec(OUT.timestamp(1)); year=ttt(1);

figure(1)
    subplot(3,2,1)
plot_Tfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
%colormap(gca,jet); 
title(['gravel road  (at year ',num2str(year),')']); grid on

figure(2)
    subplot(3,2,1)
plot_VLWCfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
title(['gravel road  (at year ',num2str(year),')']); grid on
% 
% tundra
% load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization1_output1986
% load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization1_finalState1986
% load Data\GravelRoad_PrudoeBay_realization1_output2000
% load Data\GravelRoad_PrudoeBay_realization1_finalState2000
% ttt=datevec(OUT.timestamp(1)); year=ttt(1);
% 
% figure(1)
%     subplot(3,2,2)
% plot_Tfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% colormap(gca,jet); 
% title('tundra'); grid on
% 
% figure(2)
%     subplot(3,2,2)
% plot_VLWCfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% title('tundra'); grid on
% 
% % year y
% load Data\GravelRoad_PrudoeBay_realization2_output2050
% load Data\GravelRoad_PrudoeBay_realization2_finalState2050
% ttt=datevec(OUT.timestamp(1)); year=ttt(1);
% 
% figure(1)
%     subplot(3,2,3)
% plot_Tfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% colormap(gca,jet); 
% title(['gravel road  (at year ',num2str(year),')']); grid on
% 
% figure(2)
%     subplot(3,2,3)
% plot_VLWCfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% title(['gravel road  (at year ',num2str(year),')']); grid on
% 
% tundra
% load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization1_output1986
% load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization1_finalState1986
% load Data\GravelRoad_PrudoeBay_realization1_output2050
% load Data\GravelRoad_PrudoeBay_realization1_finalState2050
% ttt=datevec(OUT.timestamp(1)); year=ttt(1);
% 
% figure(1)
%     subplot(3,2,4)
% plot_Tfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% colormap(gca,jet); 
% title('tundra'); grid on
% 
% figure(2)
%     subplot(3,2,4)
% plot_VLWCfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% title('tundra'); grid on
% % year z
% load Data\GravelRoad_PrudoeBay_realization2_output2053
% load Data\GravelRoad_PrudoeBay_realization2_finalState2053
% ttt=datevec(OUT.timestamp(1)); year=ttt(1);
% 
% figure(1)
%     subplot(3,2,5)
% plot_Tfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% colormap(gca,jet); 
% title(['gravel road  (at year ',num2str(year),')']); grid on
% 
% figure(2)
%     subplot(3,2,5)
% plot_VLWCfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% title(['gravel road  (at year ',num2str(year),')']); grid on
% 
% tundra
% load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization1_output1986
% load N:\permarisk\CryoGrid3\RunsInfrastructure\GravelRoad_PrudoeBay\GravelRoad_PrudoeBay_realization1_finalState1986
% load Data\GravelRoad_PrudoeBay_realization1_output2053
% load Data\GravelRoad_PrudoeBay_realization1_finalState2053
% ttt=datevec(OUT.timestamp(1)); year=ttt(1);
% 
% figure(1)
%     subplot(3,2,6)
% plot_Tfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% colormap(gca,jet); 
% title('tundra'); grid on
% 
% figure(2)
%     subplot(3,2,6)
% plot_VLWCfield_vs_time(OUT, FINAL.PARA, FINAL.GRID)
% title('tundra'); grid on
