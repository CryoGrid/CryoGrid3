function [T, GRID, BALANCE, TEMPORARY] = CryoGridLateralSnow( PARA, GRID, BALANCE, TEMPORARY, FORCING, T)
    
    labBarrier();
    % check preconditions
    precondition_snowExchange = checkPreconditionSnowExchange( GRID, PARA );
    if precondition_snowExchange
        fprintf('\t\t\tsync - exchanging snow\n');
        % calculate terrain index with updated surface_altitudes
        PARA.ensemble.terrain_index_snow = calculateTerrainIndexSnow(PARA.ensemble.surface_altitude, PARA.ensemble.weight);
        % calculate mobile snow
        mobile_snow = zeros( 1, numlabs );
        my_mobile_snow = 0;
        meltingConditions_index = sum(GRID.snow.Snow_w)>0 || sum(T(GRID.snow.cT_domain)>0)>0;    %snow is assumed to be immobile under melting conditions
        if ~isempty(GRID.snow.cT_domain_ub) && ~meltingConditions_index	% current realization has snow cover and no melting conditions
            i=0;
            while (abs( GRID.general.K_grid(GRID.snow.cT_domain_ub+i)-GRID.general.K_grid(GRID.snow.cT_domain_lb+1) )... 	% snow only mobile above realization-specific threshold
                    > PARA.ensemble.immobile_snow_height(labindex)+GRID.general.K_delta(GRID.snow.cT_domain_ub)) ...           % if upper cell is drifted away, immobile snow height remains
                    && (PARA.ensemble.initial_altitude(labindex)-GRID.general.K_grid(GRID.snow.cT_domain_ub+i) ...					% snow only mobile above lowermost surface altitude + snowCellSize (to prevent oscillations)
                    - (min( PARA.ensemble.surface_altitude )+GRID.snow.snowCellSize) > 1e-6 )
                my_mobile_snow = my_mobile_snow + (GRID.snow.Snow_i(GRID.snow.cT_domain_ub+i)); %only "ice" mobile
                i=i+1;
            end
        end
        mobile_snow(labindex) = my_mobile_snow;
        % exchange mobile snow amounts
        for j=1:numlabs
            if j~=labindex
                % send mobile snow amount in [m SWE]
                labSend( mobile_snow(labindex), j, 4 );
            end
        end
        for j=1:numlabs
            if j~=labindex
                % receive mobile snow amount [m SWE]
                mobile_snow(j) = labReceive(j, 4);
            end
        end
        % calculate lateral snow fluxes
        my_snow_change = calculateLateralSnowFluxes2( mobile_snow, PARA );
        % apply lateral snow fluxes directly
        if my_snow_change ~= 0
            [T, GRID] = applyLateralSnowFluxes( T, PARA, GRID, FORCING, my_snow_change );
            [GRID, T, BALANCE] = updateGRID_snow(T, GRID, PARA, BALANCE);
            BALANCE.water.dr_lateralSnow = BALANCE.water.dr_lateralSnow + my_snow_change*1000 ;
        end

        TEMPORARY.snow_flux_lateral = TEMPORARY.snow_flux_lateral + my_snow_change;
        fprintf('\t\t\tSnow flux to worker %d = %f mm SWE \n', [labindex, my_snow_change*1000] );
    end

end
