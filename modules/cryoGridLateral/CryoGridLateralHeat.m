function [ T, TEMPORARY ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid )

    labBarrier();
    % check preconditions
    precondition_heatExchange = true; %no specific conditions so far
    if precondition_heatExchange
        fprintf('\t\t\tsync - exchanging heat\n');
        % calculate lateral heat fluxes
        PACKAGE_heatExchange.T = T;
        PACKAGE_heatExchange.cT_grid = GRID.general.cT_grid;
        PACKAGE_heatExchange.k_cTgrid = k_cTgrid;
        for j=1:numlabs
            if j~=labindex
                labSend( PACKAGE_heatExchange, j, 1);
            end
        end
        % provide temporary array for lateral heat fluxes
        dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s]
        for j=1:numlabs
            if j~=labindex
                PACKAGE_heatExchange_j = labReceive(j, 1);
                [dE_dt_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/m^3/s ]
                TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + dE_dt_lateral_j .* PARA.technical.syncTimeStep.*24.*3600; % in [ J/m^3 ]
                TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  dE_dt_lateral_j .* PARA.technical.syncTimeStep.*24.*3600 .* GRID.general.K_delta ); % depth intergrated in [ J/m^2 ]
                dE_dt_lateral = dE_dt_lateral + dE_dt_lateral_j;    % sum up contributions from all realizations in [ J/m^3/s ]
            end
        end
        % apply lateral heat fluxes for entire sync interval directly
        T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600; % no division by K_delta necessary as dE_dt_lateral in [ J / m^3 / s ]
    end
end
