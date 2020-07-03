function [PARA] = updateAuxiliaryVariablesAndCommonThresholds( T, wc, GRID, PARA)    

    PARA.ensemble.surface_altitude(labindex) = getSurfaceAltitude( PARA, GRID );
    PARA.ensemble.altitude(labindex) = getAltitude( PARA, GRID );
	PARA.ensemble.water_table_altitude(labindex) = getWaterTableAltitudeFC(T, wc, GRID, PARA );% use getWaterTabelFC to account for non-saturated cells above fieldCapacity
	PARA.ensemble.soil_altitude(labindex) = getSoilAltitude( PARA, GRID );
	[ PARA.ensemble.infiltration_altitude(labindex), PARA.location.bottomBucketSoilcTIndex ] = getInfiltrationAltitude(PARA, GRID, T);

    % sending information from "labindex" to all "j"
    for j=1:numlabs
        if j~=labindex  
        %if abs(j-labindex)==1 % calc only for adjacent tiles 
        %if j-labindex==1 % calc only for adjacent tiles                        
            labSend(PARA.ensemble.surface_altitude(labindex), j, 101);
            labSend(PARA.ensemble.altitude(labindex), j, 102);
            labSend(PARA.ensemble.water_table_altitude(labindex), j, 103);
            labSend(PARA.ensemble.infiltration_altitude(labindex), j, 104);
            labSend(PARA.ensemble.soil_altitude(labindex), j, 105);
        end
    end
    % receiving ---------------------------------------------------
    % receive information from all other realizations
    for j=1:numlabs  %update surface altitudes to recalculate terrain index
        if j~=labindex
      % if abs(j-labindex)==1 % calc only for adjacent tiles
      % if j-labindex==-1 % calc only for adjacent tiles            
            PARA.ensemble.surface_altitude(j)=labReceive(j, 101);
            PARA.ensemble.altitude(j)=labReceive(j, 102);
            PARA.ensemble.water_table_altitude(j)=labReceive(j, 103);
            PARA.ensemble.infiltration_altitude(j)=labReceive(j, 104);
            PARA.ensemble.soil_altitude(j)=labReceive(j, 105);
        end
    end

    % update common thresholds for water and snow exchagne
%tsvd IS  no common threshold for snow    
%   PARA.location.absolute_maxSnow_altitude = getMaxSnowAltitude(PARA);

    %    if(PARA.IS.tile_type=='shoulder_ag')
%        % PARA.location.absolute_maxSnow_altitude = PARA.location.altitude  + 2*PARA.snow.relative_maxSnow; 
%        PARA.location.absolute_maxSnow_altitude = PARA.location.altitude  + 1.5*PARA.snow.relative_maxSnow;       
%    else
%tsvd IS   use relative maxSnow to define abs heights...! (only works if base elevation = 0m) tttt
        %PARA.location.absolute_maxSnow_altitude =  PARA.location.altitude + PARA.snow.relative_maxSnow;
        PARA.location.absolute_maxSnow_altitude =  max(0,PARA.location.altitude) + PARA.snow.relative_maxSnow; % if subsidence below 0, than define maxSnow height w.r.t. 0 baseline   ccc update if subsidence of shoulder activated...
        %PARA.location.absolute_maxSnow_altitude = PARA.snow.relative_maxSnow;
%    end
    PARA.location.absolute_maxWater_altitude = getMaxWaterAltitude(PARA);

	% update auxiliary variables in location struct
    PARA.location.altitude = PARA.ensemble.altitude(labindex);
    PARA.location.surface_altitude = PARA.ensemble.surface_altitude(labindex);
    PARA.location.soil_altitude = PARA.ensemble.soil_altitude(labindex);
    PARA.location.water_table_altitude = PARA.ensemble.water_table_altitude(labindex);
	PARA.location.infiltration_altitude = PARA.ensemble.infiltration_altitude(labindex);
    PARA.soil.infiltration_limit_altitude = PARA.location.soil_altitude - PARA.soil.infiltration_limit_depth;
