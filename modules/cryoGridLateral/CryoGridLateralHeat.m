function [ T, TEMPORARY, BALANCE ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid )
%tsvd   function [ T, TEMPORARY ] = CryoGridLateralHeat( PARA, GRID, BALANCE, TEMPORARY, T, k_cTgrid, c_cTgrid )

    dE_dt_lateral = zeros( length(GRID.general.cT_grid), 1);  %in [J/m^3/s] - temporary array for lateral heat fluxes
    F_lateral_j = zeros( length(GRID.general.cT_grid), 1); %tsvd IS
    
    labBarrier();
    % check preconditions
    precondition_heatExchange = true; %no specific conditions so far
    if precondition_heatExchange
        fprintf('\t\t\tsync - exchanging heat\n');
        % calculate lateral heat fluxes
        PACKAGE_heatExchange.T = T;
        PACKAGE_heatExchange.K_grid = GRID.general.K_grid;
        PACKAGE_heatExchange.k_cTgrid = k_cTgrid;
        PACKAGE_heatExchange.GRID.lake.cT_domain = GRID.lake.cT_domain;  %tsvd IS 

% % A) original calculation of lateral fluxes (exchange between all tiles)        
%         for j=1:numlabs
%             if j~=labindex
%                 labSend( PACKAGE_heatExchange, j, 1);
%             end
%         end       
%         for j=1:numlabs
%             if j~=labindex
%                 PACKAGE_heatExchange_j = labReceive(j, 1);
%             %    contact_length_index_j = PARA.ensemble.thermal_contact_length;  
%                 contact_length_index_j = PARA.ensemble.thermal_contact_length(j, labindex); %tsvd IS
%                 if(contact_length_index_j>0)
%                      [F_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/s ]
%         %tsvd IS     TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600; % in [ J/m^2 ]
%                      TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600; % in [ J/m^2 ]  % calculate flux through lateral cross section           
%             %tsvd    TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600 ); % depth intergrated in [ J/m^2 ]
%                     %TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600 ); % depth intergrated in [ J/m^2 ]
%                      TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + trapz(GRID.general.cT_grid, F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600 ); % depth intergrated in [ J/m^2 ]
%                     %TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(F_lateral_j ./ contact_length_index_j * PARA.technical.syncTimeStep *24.*3600 )/nansum(GRID.general.K_delta); % depth weighted sum in [ J/m^2 ]
%                %tsvd dE_dt_lateral = dE_dt_lateral + F_lateral_j./PARA.location.area;    % sum up contributions from all realizations in [ J/ m^2 /s ]
%                      dE_dt_lateral = dE_dt_lateral + F_lateral_j./ (contact_length_index_j * GRID.general.K_delta); % sum up contributions from all realizations in [ J/ m^2 /s ]
%                 end
%             end
%         end
        
% B) Adapted calculation (exchange only between adjacent tiles), fluxes calculated through vertical cross section (instead of using horizontal area), use of labSendReceive (instead of labsend and labreceive separately)
        j_prevWorker = mod(labindex - 2, numlabs) + 1; j_nextWorker = mod(labindex, numlabs) + 1; 
        j_adjacent = [j_prevWorker,j_nextWorker];
        PACKAGE_heatExchange_prev = labSendReceive(j_nextWorker,j_prevWorker,PACKAGE_heatExchange,1); 
        PACKAGE_heatExchange_next = labSendReceive(j_prevWorker,j_nextWorker,PACKAGE_heatExchange,1);
        PACKAGE_heatExchange_adjacent = [PACKAGE_heatExchange_prev,PACKAGE_heatExchange_next];
        
        for jj=1:2 % lateral heat exchange with previous and next tile
            j = j_adjacent(jj); % first exchange with previous tile, then exchange with next tile     
            PACKAGE_heatExchange_j = PACKAGE_heatExchange_adjacent(jj);
            contact_length_index_j = PARA.ensemble.thermal_contact_length(j, labindex);  
            %tsvd NOR contact_length_index_j = PARA.ensemble.thermal_contact_length(labindex);  % contact length is same for all tiles!
            
            [F_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/s ]
%             if(labindex==1 && j_prevWorker==numlabs  || labindex==numlabs && j_nextWorker==1 ) % no lateral exchange between outer tiles!
%                F_lateral_j = zeros( length(GRID.general.cT_grid), 1);
%             end
       if(contact_length_index_j >0)   
       %tsvd IS NOR     if(contact_length_index_j >0)   
%tsvd   fluxes now expressed w.r.t. vertical lateral cross section, not w.r.t horizontal "area"
     %tsvd      TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600; % in [ J/m^2 ]
            TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600; % in [ J/m^2 ]    
     %tsvd      TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600 ); % depth intergrated in [ J/m^2 ]
                %TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600 ); % depth intergrated in [ J/m^2 ]
            TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + trapz(GRID.general.cT_grid, F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600 ); % depth intergrated in [ J/m^2 ]
            %TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(F_lateral_j ./ contact_length_index_j * PARA.technical.syncTimeStep *24.*3600 )/nansum(GRID.general.K_delta); % depth weighted sum in [ J/m^2 ]
% github xice_mpi  dE_dt_lateral = dE_dt_lateral + F_lateral_j./PARA.location.area;    % sum up contributions from all realizations in [ J/ m^2 /s ]
            dE_dt_lateral = dE_dt_lateral + F_lateral_j./ (contact_length_index_j * GRID.general.K_delta); % sum up contributions from all realizations in [ J/ m^2 /s ]
%tsvd NOR   dE_dt_lateral = dE_dt_lateral + F_lateral_j; % sum up contributions from all realizations in [ J/s ]
       end
            if(numlabs==2) % if only 2 tiles, only 1 loop is needed
                break
            end
        end
% % B) Adapted calculation (exchange only between adjacent tiles), fluxes calculated through vertical cross section (instead of using horizontal area), use of labSendReceive (instead of labsend and labreceive separately)
%        j_prevWorker = mod(labindex - 2, numlabs) + 1;
%        j_nextWorker = mod(labindex, numlabs) + 1;       
%        PACKAGE_heatExchange_j = labSendReceive(j_nextWorker,j_prevWorker,PACKAGE_heatExchange,1)
%         if(labindex>1)
%             j = labindex-1; %tsvd IS  replace j by labindex-1        
% %tsvd IS    contact_length_index_j = PARA.ensemble.thermal_contact_length;  
%             contact_length_index_j = PARA.ensemble.thermal_contact_length(j, labindex);
%             [F_lateral_j, BALANCE] = calculateLateralHeatFluxes(T, k_cTgrid, PACKAGE_heatExchange_j,GRID, PARA, BALANCE, j);   % contribution from worker j to worker index in [ J/s ]
%      %tsvd      TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600; % in [ J/m^2 ]
%             TEMPORARY.dE_cell_lateral(:,j) = TEMPORARY.dE_cell_lateral(:,j) + F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600; % in [ J/m^2 ]             
%      %tsvd      TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ PARA.location.area .* PARA.technical.syncTimeStep.*24.*3600 ); % depth intergrated in [ J/m^2 ]
%                 %TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(  F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600 ); % depth intergrated in [ J/m^2 ]
%             TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + trapz(GRID.general.cT_grid, F_lateral_j ./ (contact_length_index_j * GRID.general.K_delta) * PARA.technical.syncTimeStep *24.*3600 ); % depth intergrated in [ J/m^2 ]
%                 %TEMPORARY.dE_tot_lateral(j) = TEMPORARY.dE_tot_lateral(j) + nansum(F_lateral_j ./ contact_length_index_j * PARA.technical.syncTimeStep *24.*3600 )/nansum(GRID.general.K_delta); % depth weighted sum in [ J/m^2 ]
%                 %tsvd      dE_dt_lateral = dE_dt_lateral + F_lateral_j./PARA.location.area;    % sum up contributions from all realizations in [ J/ m^2 /s ]
%             dE_dt_lateral = dE_dt_lateral + F_lateral_j./ (contact_length_index_j * GRID.general.K_delta); % sum up contributions from all realizations in [ J/ m^2 /s ]
%         end

        % apply lateral heat fluxes for entire sync interval directly
 %tsvd       T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600 ./ GRID.general.K_delta ; % division by K_delta necessary as dE_dt_lateral in [ J / m^2 / s ]
 % github xice_mpi   T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600 ./ GRID.general.K_delta ; % division by K_delta necessary as dE_dt_lateral in [ J / m^2 / s ]
       
         T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600 ./ PARA.ensemble.TileWidth(labindex); % division by tile width as dE_dt_lateral in [ J / m^2 / s ]
%tsvd IS NOR    TileVOL =  PARA.ensemble.area(labindex) * GRID.general.K_delta; 
%         if(strcmp(PARA.IS.TileType,'pile'))
%             TileVOL = TileVOL/5.; % account for different geometry for pile
%         end 
%         T = T + dE_dt_lateral./c_cTgrid.*PARA.technical.syncTimeStep.*24.*3600 ./ TileVOL; % division by grid cell volume of tile     
        % account for lateral heat fluxes in diagnostic BALANCE struct
        % (summed contribution from all connected realizations)
%tsvd  nansum is incorrect!   dE_dt_lateral is already sum over workers!       BALANCE.energy.Q_lateral = BALANCE.energy.Q_lateral + nansum(  dE_dt_lateral ) .* PARA.technical.syncTimeStep.*24.*3600;  % in [ J / m^2 ]
        BALANCE.energy.Q_lateral = BALANCE.energy.Q_lateral + dE_dt_lateral * PARA.technical.syncTimeStep.*24.*3600;  % in [ J / m^2 ]
    end
end
