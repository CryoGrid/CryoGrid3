function [precond_water] = checkPreconditionWaterExchange( T, GRID )

	precondition_waterExchange = double( T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub) ); % matches with conditions of infiltration

	for j=1:numlabs
        if j~=labindex
            labSend( precondition_waterExchange, j, 2);
        end
    end
	
    for j=1:numlabs
        if j~=labindex
            precondition_waterExchange = precondition_waterExchange + labReceive( j, 2);
        end
    end
            
	precond_water = precondition_waterExchange > 1	% at least two workers need matching conditions
