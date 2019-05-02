function fig = iPlotTemperatures(filename, OUT, PARA, GRID, cm)

    ts = OUT.timestamp();
    Ts = OUT.cryoGrid3();
    
    zs = PARA.location.initial_altitude-GRID.general.cT_grid; 
    
    fig=figure('visible','off');

    % limits
    minz = min(OUT.location.altitude - 2);
    maxz = max(OUT.location.altitude + 1);

    pcolor( ts', zs', Ts);
    hold on;
    shading flat;
    % colormap and colorbar
    caxis( [ -40, 20] );
    colormap(gca, cm.Colormap_blueautumn);
    cbar=colorbar('location','westoutside');
    % layout
    axis( [ ts(1) ts(end) minz maxz ] );
    datetick('x','mmm');%, 'keepticks'); 
    xlabel('time')
    ylabel('$z$ [m]', 'Interpreter', 'latex');
    xlabel(cbar, '$T$ [$^\circ$C]', 'Interpreter', 'latex');
    grid('on');
    
    hold off;
    
    saveas( fig, filename);
end