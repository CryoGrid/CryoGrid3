function fig = iPlotWaterContent( filename, OUT, PARA, GRID)

    ts = OUT.timestamp();
    LWCs = OUT.liquidWater();

    zs = PARA.location.initial_altitude-GRID.general.cT_grid; 

    fig=figure('visible','off');

    % limits
    minz = min(OUT.location.altitude - 2);
    maxz = max(OUT.location.altitude + 1);

    pcolor( ts', zs', LWCs);
    hold on;
    shading flat;
    % colormap and colorbar
    caxis( [ 0, 1] );
    colormap(gca, 'parula');
    colormap(gca, flipud(colormap));
    cbar=colorbar('location','eastoutside');
    % layout
    axis( [ ts(1) ts(end) minz maxz ] );
    datetick('x','mmm');%, 'keepticks'); 
    xlabel('time')
    ylabel('$z$ [m]', 'Interpreter', 'latex');
    xlabel(cbar, '$\theta_w$ [-]', 'Interpreter', 'latex');
    grid('on');

    hold off;
    
    saveas( fig, filename);
end
