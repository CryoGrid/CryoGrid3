%plot_Tfiled_tiles(numTiles,year,month,day) % year as string
% plot time slice of all tiles for indicated year, month, day  - and plot annual time series for all tiles
% first do load_TsoilData.m     (analysis_FT_WT.m for specified year)  

numTiles=5

year=2040
month=9
day=1
ti=[year month day]
    ExpLoc='Prudhoe'; 
    ExpSet='Drain10cmPond5cm'
    %ExpSet='Drain10cmPond5cm_xice'
    Scen='RCP85';
    Runs='Runs_paper/';   
load(['Data/',Runs,'/Data_',num2str(numTiles),'tiles_',num2str(year),'_',ExpSet])

%Tsoil_tiles_ti = load_TsoilData(ti);

EHbg = -1.5;
minHeight = -8; %[m]
cm = load('cm_blueautumn.mat');

%%
if(numTiles==5)
    %load(['Data/Data_FTWT_5tiles_',year])

    %checkdate=datevec(ts); idx_yearstart=find(checkdate(:,1)==year & checkdate(:,2)==1 & checkdate(:,3)==1 & checkdate(:,4)==0): idx_yearend = find(checkdate(:,1)==year+1 & checkdate(:,2)==1 & checkdate(:,3)==1 & checkdate(:,4)==0) - 1;    
%    tiles = [0 5 5 10 10 15 15 20 20 25];   
    tiles = [0 5 5 10 10 15 15 20 20 25 25 26 26 27];   
    temperature = squeeze((1:154,:)); % 
    %temperature = squeeze(Tsoil_tiles_ti(:,idx_time,:));   % Convert 3D-Space-Time-Array into 2D-Space-Array.
    temperature = [temperature(:,1), temperature(:,1), temperature(:,2), temperature(:,2), temperature(:,3), temperature(:,3),...
                   temperature(:,4), temperature(:,4), temperature(:,5), temperature(:,5),... % Double values (same values at the left and right margin on one tile).
                   NaN(154,1), NaN(154,1), temperature(:,5), temperature(:,5)]; % Double values (same values at the left and right margin on one tile). 
               z_resolution = -1 * [0:0.02:1, 1.1:0.1:10, 10.2:0.2:12.6]';   % principal height resolution  ccc subsidecne..?
    height = [z_resolution + 2.5, z_resolution + 2.5, z_resolution + 1.25, z_resolution + 1.25, z_resolution, z_resolution,...
             z_resolution, z_resolution, z_resolution, z_resolution];   % Double heights (like values) and set the absolute height depending on the specific tile.
elseif(numTiles==24)
    load(['Data/Data_FTWT_24tiles_',year])
    tiles = [0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 25 25 26];   
    temperature = squeeze(Tsoil_tiles_ti(1:154,:));   %
    temperature = [temperature(:,1), temperature(:,1), temperature(:,2), temperature(:,2), temperature(:,3), temperature(:,3),...
                   temperature(:,4), temperature(:,4), temperature(:,5), temperature(:,5), temperature(:,6), temperature(:,6),...
                   temperature(:,7), temperature(:,7), temperature(:,8), temperature(:,8), temperature(:,9), temperature(:,9),...
                   temperature(:,10), temperature(:,10), temperature(:,11), temperature(:,11), temperature(:,12), temperature(:,12),...
                   temperature(:,13), temperature(:,13), temperature(:,14), temperature(:,14), temperature(:,15), temperature(:,15),...
                   temperature(:,16), temperature(:,16), temperature(:,17), temperature(:,17), temperature(:,18), temperature(:,18),...
                   temperature(:,19), temperature(:,19), temperature(:,20), temperature(:,20), temperature(:,21), temperature(:,21),...
                   temperature(:,22), temperature(:,22), temperature(:,23), temperature(:,23), temperature(:,24), temperature(:,24),... 
                   NaN(154,1), NaN(154,1), temperature(:,24), temperature(:,24)]; % Double values (same values at the left and right margin on one tile). 
   z_resolution = -1 * [0:0.02:1, 1.1:0.1:10, 10.2:0.2:12.6]';   % principal height resolution
    height = [z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.5,...
              z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.5, z_resolution + 2.25, z_resolution + 2.25,... 
              z_resolution + 1.75, z_resolution + 1.75, z_resolution + 1.25, z_resolution + 1.25, z_resolution + 0.75, z_resolution + 0.75,... 
              z_resolution + 0.25, z_resolution + 0.25, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution,...
              z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution,...
              z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution, z_resolution,...   % Double heights (like values) and set the absolute height depending on the specific tile.
              z_resolution, z_resolution, z_resolution, z_resolution];   % Double heights (like values) and set the absolute height depending on the specific tile.
end
maxHeight = max(max(height + 1));
minTemperature = -2; maxTemperature = 1;
minTemperature = -10; maxTemperature = 5;

%%
xTicks = [0 5 10 15 20 25];
xTickLabels = {'0' '5' '10' '15' '20' '\infty'};

    figure(1) % time slice
surf(tiles, height, temperature);
view([0, 90]); axis equal; shading flat;
colormap(gca, cm.Colormap_blueautumn); cbar=colorbar('location','westoutside', 'TickDirection','out');
caxis([minTemperature, maxTemperature]); axis([tiles(1), tiles(end), minHeight, maxHeight]);
hold on
%plot([0 10],[-PARA.IS.EHbg -PARA.IS.EHbg],'--k','LineWidth',1.5) % embankment depth
%plot([10 25],[-0.3 -0.3],':k','LineWidth',1.5) % peat depth
p1 = plot3([0 10 10],[EHbg EHbg 0], [max(max(temperature)) max(max(temperature)) max(max(temperature))], '--k','LineWidth',1.5, 'DisplayName', 'Embankment'); % embankment depth
p2 = plot3([10 25],[-0.3 -0.3], [max(max(temperature)) max(max(temperature))],':k','LineWidth',1.5, 'DisplayName', 'Peat'); % peat depth
%plot3([tiles(end-1) tiles(end)],[-0.3 -0.3], [max(max(temperature)) max(max(temperature))],':k','LineWidth',1.5, 'DisplayName', 'Peat'); % peat depth
xlabel('Distance to road centre [m]'); ylabel('Depth [m]'); title(cbar, 'T[°C]');
set(gca,'box','off','color','none','TickDir','out');  
%title(datestr(time(idx_time))); hold off;
title(num2str(ti)); hold off;
legend([p1 p2], 'Location', 'northeast'); xticks(xTicks); xticklabels(xTickLabels); hold off;


%surf(time,zs,Ts)

% wenn man die Box haben will, aber die Ticks nur unten und links
%     a = gca;
%     % set box property to off and remove background color
%     set(a,'box','off','color','none')
%     % create new, empty axes with box but without ticks
%     b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
%     % set original axes as active
%     axes(a)
%     % link axes in case of zooming
%     linkaxes([a b])

% % Temperature scale stretching regarding to the whole dataset (best for comparison between timesteps).
% minTemperature = min(min([Tsoil_tiles_ti(1:116,:,1); Tsoil_tiles_ti(1:104,:,2); Tsoil_tiles_ti(1:91,:,3); Tsoil_tiles_ti(1:91,:,4); Tsoil_tiles_ti(1:91,:,5)])); 
% maxTemperature = max(max([Tsoil_tiles_ti(1:116,:,1); Tsoil_tiles_ti(1:104,:,2); Tsoil_tiles_ti(1:91,:,3); Tsoil_tiles_ti(1:91,:,4); Tsoil_tiles_ti(1:91,:,5)]));
% % Temperature scale stretching regarding to the current timestep (contrasting values within one timestep).
% minTemperature = min(min([Tsoil_tiles_ti(1:116,timestep,1); Tsoil_tiles_ti(1:104,timestep,2); Tsoil_tiles_ti(1:91,timestep,3); Tsoil_tiles_ti(1:91,timestep,4);...
%                           Tsoil_tiles_ti(1:91,timestep,5)])); 
% maxTemperature = max(max([Tsoil_tiles_ti(1:116,timestep,1); Tsoil_tiles_ti(1:104,timestep,2); Tsoil_tiles_ti(1:91,timestep,3); Tsoil_tiles_ti(1:91,timestep,4);...
%                           Tsoil_tiles_ti(1:91,timestep,5)]));
% fixed scale stretching (best for easy interpretation).
%minTemperature = -40; maxTemperature = 20;

% hotColor = [1, 0.1667, 0.1667];
% aboveZeroColor = [1, 1, 0.1667];
% belowZeroColor = [0.7829, 0.7829, 1];
% coldColor = [0.1667, 0.1667, 0.336];
% 
% if minTemperature <= 0 & maxTemperature <= 0
%     zeroRelation = -minTemperature / (maxTemperature - minTemperature);
%     customCMap = [linspace(coldColor(1), belowZeroColor(1), round(255*zeroRelation))',...
%                   linspace(coldColor(2), belowZeroColor(2), round(255*zeroRelation))',...
%                   linspace(coldColor(3), belowZeroColor(3), round(255*zeroRelation))'];
% elseif minTemperature <= 0 & maxTemperature > 0
%     zeroRelation = -minTemperature / (maxTemperature - minTemperature);
%     customCMap1 = [linspace(coldColor(1), belowZeroColor(1), round(255*zeroRelation))',...
%                    linspace(coldColor(2), belowZeroColor(2), round(255*zeroRelation))',...
%                    linspace(coldColor(3), belowZeroColor(3), round(255*zeroRelation))'];
%     customCMap2 = [linspace(aboveZeroColor(1), hotColor(1), round(255*(1 - zeroRelation)))',...
%                    linspace(aboveZeroColor(2), hotColor(2), round(255*(1 - zeroRelation)))',...
%                    linspace(aboveZeroColor(3), hotColor(3), round(255*(1 - zeroRelation)))'];
%     customCMap = [customCMap1; customCMap2];
% elseif minTemperature > 0 & maxTemperature > 0
%     disp('hot');
% end
% colormap(customCMap);
