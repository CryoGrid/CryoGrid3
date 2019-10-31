function [ T, TEMPORARY, BALANCE ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid )

    labBarrier();
    % check preconditions
    precondition_heatExchange = true; %no specific conditions so far
    if precondition_heatExchange
        fprintf('\t\t\tsync - exchanging heat\n');
        % calculate lateral heat fluxes
        PACKAGE_heatExchange.T = T;
        PACKAGE_heatExchange.K_grid = GRID.general.K_grid;
        PACKAGE_heatExchange.k_cTgrid = k_cTgrid;
        for j=1:numlabs
            if PARA.ensemble.adjacency_heat(labindex,j) % only send/recieve with connected realizations to save computation time
                labSend( PACKAGE_heatExchange, j, 10);
            end
        end
        % provide temporary array for lateral heat fluxes
        dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s]
        for j=1:numlabs
            if PARA.ensemble.adjacency_heat(labindex,j)
                PACKAGE_heatExchange_j = labReceive(j, 10);
                [F_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/s ]
                TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600; % in [ J/m^2 ]
                TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600 ); % depth intergrated in [ J/m^2 ]
                dE_dt_lateral = dE_dt_lateral + F_lateral_j./PARA.location.area;    % sum up contributions from all realizations in [ J/ m^2 /s ]
            end
        end
        % apply lateral heat fluxes for entire sync interval directly
        T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600 ./ GRID.general.K_delta ; % division by K_delta necessary as dE_dt_lateral in [ J / m^2 / s ]
        
        % account for lateral heat fluxes in diagnostic BALANCE struct
        % (summed contribution from all connected realizations)
        BALANCE.energy.Q_lateral = BALANCE.energy.Q_lateral + nansum(  dE_dt_lateral ) .* PARA.technical.syncTimeStep.*24.*3600 ; % in [ J / m^2 ]
        
    end
end
