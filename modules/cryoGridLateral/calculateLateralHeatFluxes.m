function [F, BALANCE] = calculateLateralHeatFluxes(T_index, k_index, PACKAGE_heatExchange_j, GRID, PARA, BALANCE, j)
    
    % returns a vector spanning the entire cT-grid containing the lateral heat fluxes in [J/s]
    F = zeros( length(GRID.general.cT_grid), 1);
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
        % upper and lower contact altitudes
        upper_contact_altitude = min( [ PARA.ensemble.altitude(j), PARA.ensemble.altitude(labindex) ] );
        lower_contact_altitude = max( [ altitude_Kgrid_index(end), altitude_Kgrid_j(end)] );
        % create common grid which is the union of both K-grids
        altitude_Kgrid_common = flip( union( altitude_Kgrid_index, altitude_Kgrid_j ) );
        Kdelta_common = abs( altitude_Kgrid_common(1:end-1) - altitude_Kgrid_common(2:end) );
        altitude_centerGrid_common = (altitude_Kgrid_common(1:end-1)+altitude_Kgrid_common(2:end) ) ./ 2;
        % determine the contact domain using the common centerGrid
        contactDomain = altitude_centerGrid_common < upper_contact_altitude & altitude_centerGrid_common > lower_contact_altitude;
        altitude_centerGrid_contact = altitude_centerGrid_common(contactDomain);
        
        % determine cT_indeces for common grid for both realizations (in
        % order to project back the fluxes to cells of the originial grids
        contactPointIndex_index = nan(  sum( contactDomain) ,1);
        contactPointIndex_j = nan(  sum( contactDomain), 1 );
        for k = 1:sum( contactDomain )  % there might be a moere efficient way for the following for-loop
            contactPointIndex_index(k) = find( altitude_centerGrid_contact(k) < altitude_Kgrid_index, 1, 'last' );
            contactPointIndex_j(k) = find( altitude_centerGrid_contact(k) < altitude_Kgrid_j, 1, 'last' );
        end
        
        % get T and k profiles of common grid based on original indices. no
        % interpolation is done!
        T_centerPoints_index = T_index( contactPointIndex_index );
        T_centerPoints_j = T_j( contactPointIndex_j );
        k_centerPoints_index = k_index( contactPointIndex_index );
        k_centerPoints_j = k_j( contactPointIndex_j );

        % determine effectice thermal conductivities
        k_eff = (weight_index+weight_j) ./ ( weight_index./k_centerPoints_index + weight_j./k_centerPoints_j );

        % energy flux from worker j to index in [J / s] through a lateral cross section
        F_centerPoints_index = k_eff .* (T_centerPoints_j - T_centerPoints_index) ./ distance_index_j .* contact_length_index_j .* Kdelta_common(contactDomain) ;
        
        F_accumulated_index = accumarray( contactPointIndex_index , F_centerPoints_index );
        
        F = [ F_accumulated_index; zeros( length( F ) - length( F_accumulated_index ) , 1 ) ];
              


    end
end