function [PARA] = loadSoilTypes(PARA)

	% specify one soil type per row: residualWC [%], fieldCapacity [%], alpha [1/m], n
	PARA.soil.soilTypes = [ [ 0.00, PARA.soil.fieldCapacity, 4.00, 2.0 ]; ...	% sand
							[ 0.05, PARA.soil.fieldCapacity, 0.65, 1.7 ] ];		% silt

% JAN: additional type for water? [ 0.00, 1.0, 4.00, 2.0 ]
