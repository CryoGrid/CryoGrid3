function [precond_snow] = checkPreconditionSnowExchange( GRID, PARA )

	precondition_snowExchange = double( ~isempty(GRID.snow.cT_domain_ub) && ( PARA.ensemble.surface_altitude(labindex) > PARA.ensemble.altitude(labindex)+PARA.ensemble.immobile_snow_height(labindex) ) );

	for j=1:numlabs
        if j~=labindex
            labSend( precondition_snowExchange, j, 3);
        end
    end
    for j=1:numlabs
        if j~=labindex
            precondition_snowExchange = precondition_snowExchange + labReceive( j, 3);
        end
    end

	precond_snow = precondition_snowExchange > 0;	% sufficient if one realization has mobile snow
