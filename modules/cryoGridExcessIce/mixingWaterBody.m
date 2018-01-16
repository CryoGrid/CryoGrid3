function [T] = mixingWaterBody(T, GRID)

    % mixing of temperatures in unfrozen water domain
    if GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-6 && T(GRID.soil.cT_domain_ub)>0
        unfrozenWaterBody = GRID.soil.cT_organic+GRID.soil.cT_mineral<=1e-6 & T(GRID.soil.cT_domain)>0;
        waterCells = sum( unfrozenWaterBody );
        if waterCells > 1
            Tav = sum( T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+waterCells-1) .* GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+waterCells-1) ) ...
                    ./ sum( GRID.general.K_delta(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+waterCells-1) ) ;
            T(GRID.soil.cT_domain_ub:GRID.soil.cT_domain_ub+waterCells-1) = Tav;
        end
    end

end