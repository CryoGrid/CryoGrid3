function [F, BALANCE] = calculateLateralHeatFluxes(T_index, k_index, PACKAGE_heatExchange_j, GRID, PARA, BALANCE, j)
    % returns a vector spanning the entire cT-grid containing the lateral heat fluxes in [J/s]
    %ccc index is labindex, j corresponding worker,  *_index also used for vector index...! avoid double use of infdex...
    %tsvd IS  if there is lateral heat exchange between a layer which is pond and a soil layer, the distance between these 2 tiles is re-defined (for water, i.e. T>0) assuming that the exchange is between
    %the mid-point of the soil tile and the pond border (assuming the pond being well-mixed (vertically and laterally)
    
    F = zeros( length(GRID.general.cT_grid), 1); 
    if PARA.ensemble.thermal_contact_length(labindex)>0  % calculate lateral heat flux only for laterally connected workers
        
        T_j = PACKAGE_heatExchange_j.T;
        k_j = PACKAGE_heatExchange_j.k_cTgrid;
        lakeDomain_j = PACKAGE_heatExchange_j.GRID.lake.cT_domain;  %tsvd IS

        % change to absolute altitude grid
        altitude_Kgrid_index = -GRID.general.K_grid + PARA.ensemble.initial_altitude(labindex); 
        altitude_Kgrid_j = -PACKAGE_heatExchange_j.K_grid + PARA.ensemble.initial_altitude(j);
        
        % determine the contact domain
        %tsvd IS  distance_index_j = PARA.ensemble.thermalDistance(j, labindex);
        TileWidth_index = PARA.ensemble.TileWidth(labindex) * ones(length(T_index),1);  
        TileWidth_j = PARA.ensemble.TileWidth(j) * ones(length(T_j),1);          
        TileWidth_index(GRID.lake.cT_domain & T_index>0.01) = 0.;   % cccc distance is given by distance between the centre of the adjacent soil tile and the pond border) >0.1 consider only water, not ice...  assure that no lake ice surface...!       
        TileWidth_j(lakeDomain_j & T_j>0.01) = 0.;   % distance is given by distance between the centre of the adjacent soil tile and the pond border)        
     
    switch string(PARA.Exp.Case)  %todotodo  do better in get_par_vars!?
        case 'GravelRoad'
            contact_length_index_j = PARA.ensemble.thermal_contact_length(j, labindex);  
        case 'FuelTank'
            contact_length_index_j = min(PARA.ensemble.thermal_contact_length(labindex),PARA.ensemble.thermal_contact_length(j)); 
        otherwise
            disp('no valid Experiment Case')
    end
               
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
        
        % determine cT_indeces for common grid for both realizations (in order to project back the fluxes to cells of the originial grids)
        contactPointIndex_index = nan(  sum( contactDomain) ,1);
        contactPointIndex_j = nan(  sum( contactDomain), 1 );
        for k = 1:sum( contactDomain )  % there might be a moere efficient way for the following for-loop
            contactPointIndex_index(k) = find( altitude_centerGrid_contact(k) < altitude_Kgrid_index, 1, 'last' );
            contactPointIndex_j(k) = find( altitude_centerGrid_contact(k) < altitude_Kgrid_j, 1, 'last' );
        end
        
        % get T and k profiles of common grid based on original indices. no interpolation is done!
        T_centerPoints_index = T_index( contactPointIndex_index ); 
        T_centerPoints_j = T_j( contactPointIndex_j );
        k_centerPoints_index = k_index( contactPointIndex_index );
        k_centerPoints_j = k_j( contactPointIndex_j );
        %%%weight_centerPoints_index = weight_index( contactPointIndex_index );
        %%%weight_centerPoints_j = weight_j( contactPointIndex_j );
        TileWidth_centerPoints_index = TileWidth_index( contactPointIndex_index ); 
        TileWidth_centerPoints_j = TileWidth_j( contactPointIndex_j ); 
        distance_centerPoints_index = 0.5*(TileWidth_centerPoints_index + TileWidth_centerPoints_j);
        
        k_eff = (TileWidth_centerPoints_index+TileWidth_centerPoints_j) ./ ( TileWidth_centerPoints_index./k_centerPoints_index + TileWidth_centerPoints_j./k_centerPoints_j ); %ccc
% github mpi-xice    F_centerPoints_index = k_eff .* (T_centerPoints_j - T_centerPoints_index) ./ distance_index_j .* contact_length_index_j .* Kdelta_common(contactDomain) ;
        F_centerPoints_index = k_eff .* (T_centerPoints_j - T_centerPoints_index) ./ distance_centerPoints_index .* contact_length_index_j .* Kdelta_common(contactDomain) ;   % F [W] is total flux into tile, not per unit area!   
        F_centerPoints_index(isinf(F_centerPoints_index)) = 0.; % de-activate lateral heat fluxes for adjacent grid cells which are both ponds
        F_centerPoints_index(isnan(F_centerPoints_index)) = 0.; 
        %F_accumulated_index = accumarray( contactPointIndex_index , F_centerPoints_index );       
        F_accumulated_index = accumarray( contactPointIndex_index , F_centerPoints_index);       
        F = [ F_accumulated_index; zeros( length( F ) - length( F_accumulated_index ) , 1 ) ];
        
    end
end