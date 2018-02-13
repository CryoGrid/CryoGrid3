function [PARA] = updateAuxiliaryVariablesAndCommonThresholds( T, wc, GRID, PARA, index)    

    PARA.ensemble.surface_altitude(labindex) = getSurfaceAltitude( PARA, GRID );
    PARA.ensemble.altitude(labindex) = getAltitude( PARA, GRID );
    PARA.ensemble.water_table_altitude(labindex) = getWaterTableAltitudeFC(T, wc, GRID, PARA ); %JAN: Leo uses getWaterTabelFC to account for non-saturated cells above fieldCapacity
    PARA.ensemble.active_layer_depth_altitude(labindex) = getActiveLayerDepthAltitude(PARA, GRID, T, index);
    % sending information from "labindex" to all "j"
    for j=1:numlabs
        if j~=labindex
            labSend(PARA.ensemble.surface_altitude(labindex), j, 1);
            labSend(PARA.ensemble.altitude(labindex), j, 2);
            labSend(PARA.ensemble.water_table_altitude(labindex), j, 3);
            labSend(PARA.ensemble.active_layer_depth_altitude(labindex), j, 4);
        end
    end
    % receiving ---------------------------------------------------
    % receive information from all other realizations
    for j=1:numlabs  %update surface altitudes to recalculate terrain index
        if j~=labindex
            PARA.ensemble.surface_altitude(j)=labReceive(j, 1);
            PARA.ensemble.altitude(j)=labReceive(j, 2);
            PARA.ensemble.water_table_altitude(j)=labReceive(j, 3);
            PARA.ensemble.active_layer_depth_altitude(j)=labReceive(j, 4);
        end
    end

    % update common thresholds for water and snow exchagne
    PARA.location.absolute_maxSnow_altitude = getMaxSnowAltitude(PARA);
    PARA.location.absolute_maxWater_altitude = getMaxWaterAltitude(PARA);

	% update auxiliary variables in location struct	// MORE?
    PARA.location.altitude = PARA.ensemble.altitude(labindex);
    PARA.location.surface_altitude = PARA.ensemble.surface_altitude(labindex);
    PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
	PARA.location.active_layer_depth_altitude = PARA.ensemble.active_layer_depth_altitude(labindex);
