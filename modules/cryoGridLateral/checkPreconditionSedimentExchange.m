function [precond_sediment] = checkPreconditionSedimentExchange( T, GRID, PARA )

	precondition_sedimentExchange = double( T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub) ); % matches with conditions of infiltration

	for j=1:numlabs
        if j~=labindex
            labSend( precondition_sedimentExchange, j, 41);
        end
    end
	
    for j=1:numlabs
        if j~=labindex
            precondition_sedimentExchange = precondition_sedimentExchange + labReceive( j, 41);
        end
    end
            
	precond_sediment = precondition_sedimentExchange > 1;	% at least two realizations have matching onditions
