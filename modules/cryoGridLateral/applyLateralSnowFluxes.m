function [T, GRID] = applyLateralSnowFluxes( T, PARA, GRID, FORCING, my_snow_change )
if sum(isnan(GRID.snow.Snow_i(:)))>0
   fprintf('NaN in GRID.snow.Snow_i - aLF flag1\n') 
end
if sum(isnan(GRID.snow.Snow_a(:)))>0
    fprintf('NaN in GRID.snow.Snow_a - aLF flag1\n')
end
    if ~isempty(GRID.snow.cT_domain_ub) %snow cover already exitis
        if my_snow_change>0 
            disp( 'depositing lateral snow' );
            while my_snow_change > 0
                temp_snow_flux = min( [ my_snow_change, PARA.technical.SWEperCell-GRID.snow.Snow_i(GRID.snow.cT_domain_ub) ]);

                GRID.snow.Snow_i(GRID.snow.cT_domain_ub) = GRID.snow.Snow_i(GRID.snow.cT_domain_ub) ...
                                                            + temp_snow_flux;
                if sum(isnan(GRID.snow.Snow_i(:)))>0
                    fprintf('NaN in GRID.snow.Snow_i - aLF flag2\n')
                end
                GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = GRID.snow.Snow_a(GRID.snow.cT_domain_ub) ... 
                                                            + (temp_snow_flux./(PARA.snow.rho_snow./PARA.constants.rho_w) - temp_snow_flux); % could replace this to assure correct fresh snow density
                if sum(isnan(GRID.snow.Snow_a(:)))>0
                    fprintf('NaN in GRID.snow.Snow_a - aLF flag2\n')
                end
                                                        
                % add an empty snow cell                             
                GRID.snow.cT_domain(GRID.air.cT_domain_lb)=1;
                GRID.snow.K_domain(GRID.air.K_domain_lb)=1;
                [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
                [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

                GRID.air.cT_domain(GRID.air.cT_domain_lb)=0;
                GRID.air.K_domain(GRID.air.K_domain_lb)=0; 
                [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
                [GRID.air.K_domain_lb, GRID.air.K_domain_ub] = LayerIndex(GRID.air.K_domain);

                GRID.snow.Snow_i(GRID.snow.cT_domain_ub)=0;
                if sum(isnan(GRID.snow.Snow_i(:)))>0
                    fprintf('NaN in GRID.snow.Snow_i - aLF flag3\n')
                end
                GRID.snow.Snow_w(GRID.snow.cT_domain_ub)=0;
                GRID.snow.Snow_a(GRID.snow.cT_domain_ub)=0;
                T(GRID.snow.cT_domain_ub)=T(GRID.snow.cT_domain_ub+1);
                
                if sum(isnan(GRID.snow.Snow_a(:)))>0
                    fprintf('NaN in GRID.snow.Snow_a - aLF flag3\n')
                end

                my_snow_change = my_snow_change - temp_snow_flux;

            end
        elseif my_snow_change<0
            disp( 'removing lateral snow' );
            while -my_snow_change >= GRID.snow.Snow_i(GRID.snow.cT_domain_ub) 

                my_snow_change = my_snow_change + GRID.snow.Snow_i(GRID.snow.cT_domain_ub);
                % route down water
                GRID.snow.Snow_w(GRID.snow.cT_domain_ub+1) = GRID.snow.Snow_w(GRID.snow.cT_domain_ub+1) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub); % this will not work if only one snow cell left
                % clean upper cell
                GRID.snow.Snow_i(GRID.snow.cT_domain_ub)=0;
                if sum(isnan(GRID.snow.Snow_i(:)))>0
                    fprintf('NaN in GRID.snow.Snow_i - aLF flag4\n')
                end
                GRID.snow.Snow_w(GRID.snow.cT_domain_ub)=0;
                GRID.snow.Snow_a(GRID.snow.cT_domain_ub)=0;
                T(GRID.snow.cT_domain_ub)=FORCING.i.Tair;
                if sum(isnan(GRID.snow.Snow_a(:)))>0
                    fprintf('NaN in GRID.snow.Snow_a - aLF flag4\n')
                end
                % remove upper cell
                GRID.snow.cT_domain(GRID.snow.cT_domain_ub)=0;
                GRID.snow.K_domain(GRID.snow.cT_domain_ub)=0;
                [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
                [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

                GRID.air.cT_domain(GRID.air.cT_domain_lb+1)=1;
                GRID.air.K_domain(GRID.air.K_domain_lb+1)=1;
                [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
                [GRID.air.K_domain_lb, GRID.air.K_domain_ub] =   LayerIndex(GRID.air.K_domain);

            end

            % remove remaining snow_flux from upper cell
            GRID.snow.Snow_i(GRID.snow.cT_domain_ub) = GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + my_snow_change;
            if sum(isnan(GRID.snow.Snow_i(:)))>0
                fprintf('NaN in GRID.snow.Snow_i - aLF flag5\n')
            end
            %GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = GRID.snow.Snow_a(GRID.snow.cT_domain_ub) + (my_snow_change./(PARA.snow.rho_snow./PARA.constants.rho_w) - my_snow_change);
            if sum(isnan(GRID.snow.Snow_a(:)))>0
                fprintf('NaN in GRID.snow.Snow_a - aLF flag5\n')
            end
            GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = GRID.snow.Snow_i(GRID.snow.cT_domain_ub) .* ( (PARA.constants.rho_w ./ PARA.snow.rho_snow ) - 1 ); % this assures fresh snow density for upper cell
            if sum(isnan(GRID.snow.Snow_a(:)))>0
                fprintf('NaN in GRID.snow.Snow_a - aLF flag6\n')
            end
        else % lateral_snow_flux==0
            %do nothing
        end

    elseif my_snow_change>0 %no snow cover and positive lateral input
        %---------- add the new snow into initial SWE variable in case of no snow cover------------------
        GRID.snow.SWEinitial = GRID.snow.SWEinitial + my_snow_change;
    elseif my_snow_change<0 % no snow cover, but negative lateral flow calculated
        warning('snow exchange - negative lateral snow flux without existing snow cover');
    end




end