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
%tsvd allow for heat exchange over full vertical profile    contact_altitude = min( [ PARA.ensemble.altitude(j), PARA.ensemble.altitude(index) ] ) - distance_index_j;  %below this depth, grid cells will exchange heat
        contact_altitude = min( [ PARA.ensemble.altitude(j), PARA.ensemble.altitude(index) ] );  %below this depth, grid cells will exchange heat
        contact_domain = altitude_cTgrid_index <= contact_altitude;  %all cells in the current ensemble member

        % interpolate j-values to index-grid   
        altitude_cTgrid_j = -PACKAGE_heatExchange_j.cT_grid + PARA.ensemble.initial_altitude(j);
        min_contact_altitude = min( [altitude_cTgrid_index(end), altitude_cTgrid_j(end)] ); %zzz max instead of min   def. also max_contact_alt ?
        altitude_cTgrid_index(end) = min_contact_altitude;
        altitude_cTgrid_j(end) = min_contact_altitude;
        T_j = PACKAGE_heatExchange_j.T;
        k_j = PACKAGE_heatExchange_j.k_cTgrid;
        if ~isreal(altitude_cTgrid_index); disp('altitude_cTgrid_index contains complex values');  end  %likley not relevant anymore
        if ~isreal(altitude_cTgrid_j);     disp('altitude_cTgrid_j contains complex values');      end
        if ~isreal(T_j);                   disp('T_j contains complex values');                    end
        if ~isreal(k_j);                   disp('k_j contains complex values');                    end
        
        T_interp_j = interp1( altitude_cTgrid_j, T_j, altitude_cTgrid_index, 'linear'); 
%         try
%             assert( sum( isnan( T_interp_j ) )==0, 'calc lat heat fluxes - error in T interpolation') %ttt
%         catch
%            save Data_Tinterp T_interp_j altitude_cTgrid_j T_j altitude_cTgrid_index index 
% %           error('interpolation NAN T error')
%         end
        k_interp_j = interp1( altitude_cTgrid_j, k_j, altitude_cTgrid_index, 'linear');
%         try
%             assert( sum( isnan( k_interp_j ) )==0, 'calc lat heat fluxes - error in k interpolation') %ttt
%         catch
%             save Data_kinterp k_interp_j altitude_cTgrid_j k_j altitude_cTgrid_index index
%             error('interpolation NAN K error')
%         end
        % determine effectice thermal conductivities
        k_eff = (weight_index+weight_j) ./ ( weight_index./k_index + weight_j./k_interp_j );

        % energy flux from worker j to index in [J / m^3 / s]
        dE_dt(contact_domain) = k_eff(contact_domain) .* (T_interp_j(contact_domain)-T_index(contact_domain)) ./ distance_index_j .* contact_length_index_j ./ PARA.ensemble.area(index) ; % in [ J / m^3 / s ]
              
%        T_index(contact_domain) = T_index(contact_domain) +  dE_dt_j ./ c_index(contact_domain) .* PARA.technical.syncTimeStep .* 24 .* 3600;
 
        % balance is not correct
       %tsvd  BALANCE.energy.Q_lateral(contact_domain) = BALANCE.energy.Q_lateral(contact_domain) + dE_dt(contact_domain) .* PARA.technical.syncTimeStep .* 24 .* 3600 .* GRID.general.K_delta(contact_domain); % in [ J / m^2 ]
        BALANCE.energy.Q_lateral(contact_domain) = k_eff(contact_domain) .* (T_interp_j(contact_domain)-T_index(contact_domain)) ./ distance_index_j; % in W/m2 
        
    end
end
