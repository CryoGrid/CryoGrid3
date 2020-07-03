function plot_FT_WT_tiles(year)
% plot frost and water tables for all tiles
% first generate data (using plot_FT_WT.m) or load Data

% 2 tiles
load(['Data/Data_FTWT_2tiles_',year]) 
distanceRC(2)=24; tileBoarders(2)=24;
distanceRC_2tiles=distanceRC; tileBoarders_2tiles=tileBoarders;
FT_2tiles=FT_tiles; WT_2tiles=WT_tiles;

% 5 tiles
load(['Data/Data_FTWT_5tiles_',year]) 
%disp('reset last tile width....!!!')
%distanceRC(5)=22.5; tileBoarders(5)=25
distanceRC_5tiles=distanceRC; tileBoarders_5tiles=tileBoarders;
FT_5tiles=FT_tiles; WT_5tiles=WT_tiles;

% 24 tiles
load(['Data/Data_FTWT_24tiles_',year]) 
%load(['Data_FTWT_5tiles_equiD_',year]) % 5 tiles equidistant
distanceRC_24tiles=distanceRC; tileBoarders_24tiles=tileBoarders;
FT_24tiles=FT_tiles; WT_24tiles=WT_tiles;

if(str2num(year)>2074); zmin=-10; elseif(str2num(year)>2049) ; zmin=-8; else; zmin=-6; end

figure(33)
    subplot(2,2,1)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Mineral soil');
hold on
plot(distanceRC_2tiles,FT_2tiles(ind_Jun1,:),'xr','DisplayName','Frost table')
plot(distanceRC_2tiles,WT_2tiles(ind_Jun1,:),'xb','DisplayName','Water table') %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
plot(distanceRC_5tiles,FT_5tiles(ind_Jun1,:),'*r','DisplayName','Frost table','LineWidth',2)
plot(distanceRC_5tiles,WT_5tiles(ind_Jun1,:),'*b','DisplayName','Water table','LineWidth',2) %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
plot(distanceRC_24tiles,FT_24tiles(ind_Jun1,:),'+-r','DisplayName','Frost table')
plot(distanceRC_24tiles,WT_24tiles(ind_Jun1,:),'+-b','DisplayName','Water table') %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
plot(distanceRC,nanmax(snowHeight_tiles)+Soil_surf_tiles(1,:),'+-c','DisplayName','Max Snow height') % maximum snow height in year
grid on; axis([0 tileBoarders_24tiles(end) zmin 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of June ',year]); 
xxx=xlim; yyy=ylim; text(xxx(1),yyy(2),[ExpLoc,' ',ExpSet],'fontsize',8,'Interpreter','none')
%legend('show', 'Location', 'southeast'); % legend('frost table','water table','soil surface')
    subplot(2,2,2)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Mineral soil');
hold on
plot(distanceRC_2tiles,FT_2tiles(ind_Jul1,:),'xr','DisplayName','Frost table')
plot(distanceRC_2tiles,WT_2tiles(ind_Jul1,:),'xb','DisplayName','Water table')
plot(distanceRC_5tiles,FT_5tiles(ind_Jul1,:),'*r','DisplayName','Frost table','LineWidth',2)
plot(distanceRC_5tiles,WT_5tiles(ind_Jul1,:),'*b','DisplayName','Water table','LineWidth',2) %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
plot(distanceRC_24tiles,FT_24tiles(ind_Jul1,:),'+-r','DisplayName','Frost table')
plot(distanceRC_24tiles,WT_24tiles(ind_Jul1,:),'+-b','DisplayName','Water table') %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
grid on; axis([0 tileBoarders_24tiles(end) zmin 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of July ',year]); 
xxx=xlim; yyy=ylim; text(xxx(1),yyy(2),[ExpLoc,' ',ExpSet],'fontsize',8,'Interpreter','none')
%legend('show', 'Location', 'southeast'); 
    subplot(2,2,3)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Mineral soil');
hold on
plot(distanceRC_2tiles,FT_2tiles(ind_Aug1,:),'xr','DisplayName','Frost table')
plot(distanceRC_2tiles,WT_2tiles(ind_Aug1,:),'xb','DisplayName','Water table')
plot(distanceRC_5tiles,FT_5tiles(ind_Aug1,:),'*r','DisplayName','Frost table','LineWidth',2)
plot(distanceRC_5tiles,WT_5tiles(ind_Aug1,:),'*b','DisplayName','Water table','LineWidth',2) %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
plot(distanceRC_24tiles,FT_24tiles(ind_Aug1,:),'+-r','DisplayName','Frost table')
plot(distanceRC_24tiles,WT_24tiles(ind_Aug1,:),'+-b','DisplayName','Water table') %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
grid on; axis([0 tileBoarders_24tiles(end) zmin 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of Aug ',year]); 
xxx=xlim; yyy=ylim; text(xxx(1),yyy(2),[ExpLoc,' ',ExpSet],'fontsize',8,'Interpreter','none')
%legend('show', 'Location', 'southeast'); 
    subplot(2,2,4)
patch('Faces', snowFaces, 'Vertices', snow, 'FaceColor', snowColors, 'FaceAlpha', 0.7, 'DisplayName', 'Snow_{Max}'); patch('Faces', peatFaces, 'Vertices', peat, 'FaceColor', peatColors, 'FaceAlpha', 0.7, 'DisplayName', 'Peat');
patch('Faces', gravelFaces, 'Vertices', gravel, 'FaceColor', gravelColors, 'FaceAlpha', 0.7, 'DisplayName', 'Gravel'); patch('Faces', mineralsoilFaces, 'Vertices', mineralsoil, 'FaceColor', mineralsoilColors, 'FaceAlpha', 0.7, 'DisplayName', 'Mineral soil');
hold on
plot(distanceRC_2tiles,FT_2tiles(ind_Sep1,:),'xr','DisplayName','Frost table')
plot(distanceRC_2tiles,WT_2tiles(ind_Sep1,:),'xb','DisplayName','Water table')
plot(distanceRC_5tiles,FT_5tiles(ind_Sep1,:),'*r','DisplayName','Frost table','LineWidth',2)
plot(distanceRC_5tiles,WT_5tiles(ind_Sep1,:),'*b','DisplayName','Water table','LineWidth',2) %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
plot(distanceRC_24tiles,FT_24tiles(ind_Sep1,:),'+-r','DisplayName','Frost table')
plot(distanceRC_24tiles,WT_24tiles(ind_Sep1,:),'+-b','DisplayName','Water table') %,x_GRprofile,y_GRprofile,'+--k',[0 15],[-1.5 -1.5],'--k')
grid on; axis([0 tileBoarders_24tiles(end) zmin 3]); 
xlabel('Distance from Road Centre [m]'); ylabel('frost and water table [m]'); title(['1st of Sep ',year]); 
xxx=xlim; yyy=ylim; text(xxx(1),yyy(2),[ExpLoc,' ',ExpSet],'fontsize',8,'Interpreter','none')
%legend('show', 'Location', 'southeast'); 

set(gcf,'units','normalized','outerposition',[0 0 1 1])
%set(figure, 'position', [0, 0, 800, 200]) 
saveas(gcf,['FTWTtiles_',year],'tif')
