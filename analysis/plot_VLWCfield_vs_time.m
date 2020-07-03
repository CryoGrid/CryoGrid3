function fig = plot_VLWCfield_vs_time(OUT, PARA, GRID)

    ts = OUT.timestamp();
    LWCs = OUT.liquidWater();
    zs = PARA.location.initial_altitude-GRID.general.cT_grid; 
    %fig=figure('visible','off');

    % limits
    minz = min(OUT.location.altitude - 2); maxz = max(OUT.location.altitude + 1);

    pcolor( ts', zs', LWCs);
    
    hold on; shading flat;
    plot([ts(1) ts(end)],[0 0],'k:',[ts(1) ts(end)],[-0.5 -0.5],'k--',[ts(1) ts(end)],[-0.8 -0.8],'k')

    % colormap and colorbar
    caxis( [ 0, 1] ); colormap(gca, 'parula'); colormap(gca, flipud(colormap)); cbar=colorbar('location','westoutside');
    % layout
        axis( [ ts(1) ts(end) -1.5 2 ] ); %axis( [ ts(1) ts(end) minz maxz ] ); 
    datetick('x','mmm');%, 'keepticks'); 
    ylabel('z [m]'); title(cbar, 'LWC'); grid on
    hold off;