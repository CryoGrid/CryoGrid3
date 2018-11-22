function [dE_dt, BALANCE] = calculateLateralHeatFluxes(T_index, k_index, PACKAGE_heatExchange_j, GRID, PARA, BALANCE, j)
    
    dE_dt = zeros( length(GRID.general.cT_grid), 1);
    if PARA.ensemble.thermal_contact_length(labindex,j)>0  % calculate lateral heat flux only for laterally connected workers
        
        T_j = PACKAGE_heatExchange_j.T;
        k_j = PACKAGE_heatExchange_j.k_cTgrid;
        
        % change to absolute altitude grid
        altitude_Kgrid_index = -GRID.general.K_grid + PARA.ensemble.initial_altitude(labindex);
        altitude_Kgrid_j = -PACKAGE_heatExchange_j.K_grid + PARA.ensemble.initial_altitude(j);
        
        
        % determine the contact domain
        distance_index_j = PARA.ensemble.thermalDistance(j, labindex);
        weight_index = PARA.ensemble.weight(labindex);
        weight_j = PARA.ensemble.weight(j);
        contact_length_index_j = PARA.ensemble.thermal_contact_length(j, labindex);
        
        %contact_altitude = min( [ PARA.ensemble.altitude(j), PARA.ensemble.altitude(labindex) ] );  %below this depth, grid cells will exchange heat
        %contact_domain = altitude_cTgrid_index <= contact_altitude;  %all cells in the current ensemble member
        
        

        upper_contact_altitude = min( [ PARA.ensemble.altitude(j), PARA.ensemble.altitude(labindex) ] );
        lower_contact_altitude = max( [ altitude_Kgrid_index(end), altitude_Kgrid_j(end)] );
        
        altitude_Kgrid_common = flip( union( altitude_Kgrid_index, altitude_Kgrid_j ) );   
        altitude_centerGrid_common = (altitude_Kgrid_common(1:end-1)+altitude_Kgrid_common(2:end) ) ./ 2;
        
        contactDomain_centerGrid_common = altitude_centerGrid_common < upper_contact_altitude & altitude_centerGrid_common > lower_contact_altitude;

        % determine cT_index for common grid for both realizations
        centerPointIndex_common_index = nan(  length( altitude_centerGrid_common ) ,1);
        centerPointIndex_common_j = nan(  length( altitude_centerGrid_common ), 1 );
        % this is super inefficient and should be replaced by some smart vector operation
        for k = 1:length(altitude_centerGrid_common)
            centerPointIndex_common_index(k) = find( altitude_centerGrid_common < altitude_Kgrid_index, 1, 'last' );
            centerPointIndex_common_j(k) = find( altitude_centerGrid_common < altitude_Kgrid_j, 1, 'last' );
        end
        
        T_centerPoints_index = T_index( centerPointIndex_common_index(contactDomain_centerGrid_common ) );
        T_centerPoints_j = T_j( centerPointIndex_common_j(contactDomain_centerGrid_common) );
        
        k_centerPoints_index = k_index( centerPointIndex_common_index(contactDomain_centerGrid_common ) );
        k_centerPoints_j = k_j( centerPointIndex_common_j(contactDomain_centerGrid_common ) );
        

        % determine effectice thermal conductivities
        k_eff = (weight_index+weight_j) ./ ( weight_index./k_centerPoints_index + weight_j./k_centerPoints_j );

        % energy flux from worker j to index in [J / m^3 / s]
        dE_dt_centerPoints_index = k_eff .* (T_centerPoints_j - T_centerPoints_index) ./ distance_index_j .* contact_length_index_j ./ PARA.ensemble.area(labindex) ;
        
        dE_dt_accumulated_index = accumarray( centerPointIndex_common_index(contactDomain_centerGrid_common), dE_dt_centerPoints_index );
        
        dE_dt = [ dE_dt_accumulated_index; zeros( length( dE_dt ) - length( dE_dt_accumulated_index ) , 1 ) ];
              

        % balance is not correct
        %BALANCE.energy.Q_lateral(contact_domain) = BALANCE.energy.Q_lateral(contact_domain) + dE_dt(contact_domain) .* PARA.technical.syncTimeStep .* 24 .* 3600 .* GRID.general.K_delta(contact_domain); % in [ J / m^2 ]
        
    end
end