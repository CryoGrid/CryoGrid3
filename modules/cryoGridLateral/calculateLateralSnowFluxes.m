function snow_flux_j = calculateLateralSnowFluxes( mobile_snow_index, PACKAGE_snowExchange_j, GRID, PARA, index, j)

    snow_flux_j = 0;
    snow_diffusivity = PARA.ensemble.snow_diffusivity; % in [m^2/s] constant of proportionality
    
    
    distance_index_j = PARA.ensemble.distanceBetweenPoints(index,j);
    surface_altitude_index = PARA.ensemble.surface_altitude(index);
    surface_altitude_j = PARA.ensemble.surface_altitude(j);
    snow_w_j = PACKAGE_snowExchange_j.snow.Snow_w;
    mobile_snow_j = PACKAGE_snowExchange_j.mobile_snow;
    
    % preconditions: points connected and minimum surface elevation difference in order to prevent oscillations between two workers  
    if distance_index_j~=0 && abs(surface_altitude_index-surface_altitude_j) > GRID.snow.snowCellSize         
    
            
            hasMobileSnow_index = surface_altitude_index > PARA.ensemble.altitude(index)+PARA.ensemble.immobile_snow_height(index);
            hasMobileSnow_j     = surface_altitude_j > PARA.ensemble.altitude(j)+PARA.ensemble.immobile_snow_height(j);
            
            isMelting_index = sum(GRID.snow.Snow_w)>0;
            isMelting_j = sum(snow_w_j) > 0;
            
            if isMelting_index || isMelting_j
                disp('melting conditions');
            end
            
            % this does not account for wind speed so far
            maxSnowFlux = snow_diffusivity .* abs(surface_altitude_j - surface_altitude_index) ./ distance_index_j .* PARA.ensemble.snow_contact_length(index,j) ./ PARA.ensemble.area(index); % in [m SWE / s]
       
            
               
            
            if surface_altitude_j > surface_altitude_index && hasMobileSnow_j && ~isMelting_j  % worker index receiving drift snow from j
                
                snow_flux_j = min( [maxSnowFlux, mobile_snow_j * PARA.ensemble.weight(j) / PARA.ensemble.weight(index) ]);
                
            elseif surface_altitude_index > surface_altitude_j && hasMobileSnow_index && ~isMelting_index   % worker index is depositing drift snow to j
                
                snow_flux_j = -min( [ maxSnowFlux, mobile_snow_index ] ) ;
                
            else % no drift snow exchange
                
                snow_flux_j = 0;
                
            end
                    
    end
    
end