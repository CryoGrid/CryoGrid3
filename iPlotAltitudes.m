function [ ] = iPlotAltitudes( filename, OUT, PARA )

    ts = OUT.timestamp();
    water_table = OUT.location.water_table_altitude();
    frost_table = OUT.location.infiltration_altitude();
    surface_level = OUT.location.surface_altitude;
    terrain_level = OUT.location.altitude();
    try
        soil_level = OUT.location.soil_altitude();
    catch
        soil_level = nanmin( OUT.soil.topPosition, OUT.soil.lakeFloor) + PARA.location.initial_altitude;
    end

    fig=figure('visible','off');
    
    p_area = area( ts, [ soil_level, terrain_level-soil_level, surface_level-terrain_level] );
    hold on;
    p_area(1).FaceColor = [210 180 140]./255;
    p_area(1).FaceAlpha = 0.5;
    p_area(1).DisplayName = 'soil domain';
    p_area(2).FaceColor = [176 196 222]./255;
    p_area(2).FaceAlpha = 0.5;
    p_area(2).DisplayName = 'lake domain';
    p_area(3).FaceColor = [220 220 220]./255;
    p_area(3).FaceAlpha = 0.5;
    p_area(3).DisplayName = 'snow domain';
    
    
    plot(ts, water_table, 'blue', 'LineWidth', 2, 'DisplayName','water table');
    plot(ts, frost_table, 'red', 'LineWidth', 2, 'DisplayName','frost table');

    axis( [ ts(1) ts(end) min(soil_level)-1.5 max(surface_level)+0.5 ] );
    %ax=gca;
    %ax.XTick = linspace(ts(1), ts(end), 13);
    datetick('x','mmm');     
    xlabel('time');
    ylabel('altitude [m asl]');
	title(datestr(ts(1), 'yyyy'));
    legend('show', 'Location', 'southeast');
    grid('on');
    
    hold off;
    
    saveas( fig, filename);

end

