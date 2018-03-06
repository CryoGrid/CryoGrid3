function [PARA] = loadSoilTypes( PARA )

	% specify one soil type per row: residualWC [%], fieldCapacity [%], alpha [1/m], n
	PARA.soil.soilTypes = [ [ 0.00, PARA.soil.fieldCapacity, 4.00, 2.0 ]; ...	% sand
							[ 0.05, PARA.soil.fieldCapacity, 0.65, 1.7 ]; ...	% silt
                            [ 0.00, 0.00                   , 4.00, 2.0 ] ];     % pond (freeze curve of sand but fieldCapacity = 0)