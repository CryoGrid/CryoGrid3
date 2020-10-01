function [PARA] = loadSoilTypes( PARA )
% specify one soil type per row: residualWC [%], fieldCapacity [%], alpha [1/m], n (determine freeze curve)
%tsvd  residual WC set to 5% for sand, peat, and gravel

    PARA.soil.soilTypes = [ [ 0.05, PARA.soil.fieldCapacity_MS, 4.00, 2.0 ]; ...	% 1: sand
							[ 0.05, PARA.soil.fieldCapacity_MS, 0.65, 1.7 ]; ...	% 2: silt
                            [ 0.00, 0.00                   , 4.00, 2.0 ]; ...       % 3: pond (freeze curve of sand but fieldCapacity = 0)
%tsvd IS introduce new soil types 'Peat' and 'Gravel' with specified field capacities and freeze curve parameters for sand
                            [ 0.05, PARA.soil.fieldCapacity_Peat, 4.00, 2.0 ]; ...            % 4: Peat 
                            [ 0.05, PARA.soil.fieldCapacity_Gravel, 4.00, 2.0 ]; ...	      % 5: Gravel                       
                            [ 0.05, PARA.soil.fieldCapacity_Gravel_surface, 4.00, 2.0 ]; ...  % 6: surface Gravel                       
                            [ 0.00, 0.01, 4.00, 2.0 ]; ...                                    % 7: concrete     %tsvd NOR                  
                            [ 0.00, 0.01, 4.00, 2.0 ]; ...                                    % 8: steel                 
                            [ 0.00, 0.01, 4.00, 2.0 ] ];                                      % 9: styrofoam    so far dummy values...      