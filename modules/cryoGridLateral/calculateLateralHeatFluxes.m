function [dE_dt, BALANCE] = calculateLateralHeatFluxes(T_index, k_index, PACKAGE_heatExchange_j, GRID, PARA, BALANCE, j)
    
	index = labindex;
    dE_dt = zeros( length(GRID.general.cT_grid), 1);
    if PARA.ensemble.thermal_contact_length(index,j)>0  % calculate lateral heat flux only for laterally connected workers
        
        % change to absolute altitude grid
        altitude_cTgrid_index = -GRID.general.cT_grid + PARA.ensemble.initial_altitude(index);

        % determine the contact domain
        distance_index_j = PARA.ensemble.distanceBetweenPoints(j, index);
        weight_index = PARA.ensemble.weight(index);
        weight_j = PARA.ensemble.weight(j);
        contact_length_index_j = PARA.ensemble.thermal_contact_length(j, index);
        contact_altitude = min( [ PARA.ensemble.altitude(j), PARA.ensemble.altitude(index) ] ) - distance_index_j;  %below this depth, grid cells will exchange heat
        contact_domain = altitude_cTgrid_index <= contact_altitude;  %all cells in the current ensemble member

        % interpolate j-values to index-grid   
        altitude_cTgrid_j = -PACKAGE_heatExchange_j.cT_grid + PARA.ensemble.initial_altitude(j);
        min_contact_altitude = min( [altitude_cTgrid_index(end), altitude_cTgrid_j(end)] );
        altitude_cTgrid_index(end) = min_contact_altitude;
        altitude_cTgrid_j(end) = min_contact_altitude;
        T_j = PACKAGE_heatExchange_j.T;
        k_j = PACKAGE_heatExchange_j.k_cTgrid;
        if ~isreal(altitude_cTgrid_index) %likley not relevant anymore
            disp('altitude_cTgrid_index contains complex values');
        end
        if ~isreal(altitude_cTgrid_j)
            disp('altitude_cTgrid_j contains complex values');
        end
        if ~isreal(T_j)
            disp('T_j contains complex values');
        end
        if ~isreal(k_j)
            disp('k_j contains complex values');
        end
        T_interp_j = interp1( altitude_cTgrid_j, T_j, altitude_cTgrid_index, 'linear');
        k_interp_j = interp1( altitude_cTgrid_j, k_j, altitude_cTgrid_index, 'linear');

        % determine effectice thermal conductivities
        k_eff = (weight_index+weight_j) ./ ( weight_index./k_index + weight_j./k_interp_j );

        % energy flux from worker j to index in [J / m^3 / s]
        dE_dt(contact_domain) = k_eff(contact_domain) .* (T_interp_j(contact_domain)-T_index(contact_domain)) ./ distance_index_j .* contact_length_index_j ./ PARA.ensemble.area(index) ; % in [ J / m^3 / s ]
              
%        T_index(contact_domain) = T_index(contact_domain) +  dE_dt_j ./ c_index(contact_domain) .* PARA.technical.syncTimeStep .* 24 .* 3600;

        % balance is not correct
        BALANCE.energy.Q_lateral(contact_domain) = BALANCE.energy.Q_lateral(contact_domain) + dE_dt(contact_domain) .* PARA.technical.syncTimeStep .* 24 .* 3600 .* GRID.general.K_delta(contact_domain); % in [ J / m^2 ]
        
    end
end
