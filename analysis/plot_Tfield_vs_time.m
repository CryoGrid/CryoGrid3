%function fig = plot_Tfield_vs_time(OUT, PARA, GRID, cm)
function fig = plot_Tfield_vs_time(OUT, PARA, GRID)

load cm_blueautumn.mat

    ts = OUT.timestamp();
%    Ts = OUT.cryoGrid3();
%    zs = PARA.location.initial_altitude-GRID.general.cT_grid; 
    Ts = squeeze(Tsoil_tiles(:,:,1));
    zs = zsoil_tiles(:,1);
    
%    fig=figure('visible','off');
    % limits
    minz = min(OUT.location.altitude - 6.5); maxz = max(OUT.location.altitude + 1);

    pcolor( ts', zs', Ts);    
    hold on; shading flat;
    plot([ts(1) ts(end)],[0 0],'k:',[ts(1) ts(end)],[-0.5 -0.5],'k--',[ts(1) ts(end)],[-0.8 -0.8],'k')
%    caxis( [ -40, 20] );  %caxis( [ -30, 20] );  shift 0° border....! 
    caxis( [ -20, 10] );  %caxis( [ -30, 20] );  shift 0° border....! 

%    colormap(gca, cm.Colormap_blueautumn); 
    colormap(gca, Colormap_blueautumn);  cbar=colorbar('location','westoutside');
%    axis( [ ts(1) ts(end) minz maxz ] ); 
    axis( [ ts(1) ts(end) -1.5 2 ] ); 
    datetick('x','mmm');%, 'keepticks'); 
    ylabel('z [m]');   %ylabel('$z$ [m]', 'Interpreter', 'latex');%    xlabel(cbar, '$T$ [$^\circ$C]', 'Interpreter', 'latex');
    title(cbar, 'T [°C]'); grid on 
    hold off;

%     figure(3)
%     h=imagesc(ts,zs,Ts);
%     set(h,'AlphaData',~isnan(Ts)); set(gca,'Ydir','normal')
