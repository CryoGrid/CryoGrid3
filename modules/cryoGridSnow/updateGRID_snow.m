function [GRID, T, BALANCE] = updateGRID_snow(T, GRID, PARA, BALANCE)

    snowCellSize=GRID.snow.snowCellSize;

    if isempty(GRID.snow.cT_domain_lb)==1 %no snow exists

        if GRID.snow.SWEinitial>=PARA.technical.SWEperCell/2   %create and initialize first cell of snow

            %------ modify snow and air grid -----------------------------
            GRID.snow.cT_domain(GRID.air.cT_domain_lb)=1;
            GRID.snow.K_domain(GRID.air.K_domain_lb)=1;
            [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
            [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

            GRID.air.cT_domain(GRID.air.cT_domain_lb)=0;
            GRID.air.K_domain(GRID.air.K_domain_lb)=0; 
            [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
            [GRID.air.K_domain_lb, GRID.air.K_domain_ub] =   LayerIndex(GRID.air.K_domain);

            % ------- update SWE grid -------------------------------------
            GRID.snow.Snow_i(GRID.snow.cT_domain_ub) = GRID.snow.SWEinitial;
            GRID.snow.Snow_w(GRID.snow.cT_domain_ub) = 0;
            GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = (GRID.snow.SWEinitial./(PARA.snow.rho_snow./PARA.constants.rho_w) - GRID.snow.SWEinitial);      
            GRID.snow.SWEinitial=0;       

            % -------- update K grid -------------------------------------
            assert(~isnan(GRID.general.K_grid(GRID.snow.cT_domain_ub+1)),' error in uppermost lake cell K_grid'); %tsvd
            
            GRID.general.K_grid(GRID.snow.cT_domain_ub) = GRID.general.K_grid(GRID.snow.cT_domain_ub+1)  - ( GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + GRID.snow.Snow_a(GRID.snow.cT_domain_ub) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub) );
    %tsvd        GRID.general.K_grid(GRID.snow.K_domain_ub)=-1.*( GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + GRID.snow.Snow_a(GRID.snow.cT_domain_ub) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub) );
            T(GRID.snow.cT_domain_ub)=T(GRID.air.cT_domain_lb);
            assert( sum(isnan( GRID.general.K_grid ) )==0, 'updateGRID_snow - error xxx')
            assert( sum(isnan(GRID.snow.Snow_i) )==0,' GRID.snow.snow_i NAN  1') 
        end

    else   %snow exists

        check_change=false;  
        assert( sum(isnan(GRID.snow.Snow_i) )==0,' GRID.snow.snow_i NAN  1') 
       
        GRID.general.K_grid(GRID.snow.cT_domain_ub) = GRID.general.K_grid(GRID.snow.cT_domain_ub+1) -...
            ( GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub) + GRID.snow.Snow_a(GRID.snow.cT_domain_ub)); %updates the position of the uppermost snow grid cell                
            assert( ~isnan( GRID.general.K_grid(GRID.snow.cT_domain_ub) ), 'updateGRID_snow - error in uppermost snow cell position' );
            if(~isempty(GRID.lake.water.cT_domain_ub)) 
                assert(GRID.lake.water.cT_domain(GRID.lake.water.cT_domain_ub)+GRID.snow.cT_domain(GRID.lake.water.cT_domain_ub-1)<2,'snow on lake!');
                % assert(GRID.lake.water.cT_domain_ub-GRID.snow.cT_domain_lb>1,'snow on lake in updateGRID_snow.m!')
            end
        
        if GRID.snow.Snow_i(GRID.snow.cT_domain_ub)>=1.5.*PARA.technical.SWEperCell  %create new grid cell

            %------ modify snow and air grid -----------------------------
            GRID.snow.cT_domain(GRID.air.cT_domain_lb)=1;
            GRID.snow.K_domain(GRID.air.K_domain_lb)=1;
            [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
            [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

            GRID.air.cT_domain(GRID.air.cT_domain_lb)=0;
            GRID.air.K_domain(GRID.air.K_domain_lb)=0; 
            [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
            [GRID.air.K_domain_lb, GRID.air.K_domain_ub] =   LayerIndex(GRID.air.K_domain);

            % ------- update SWE grid
            GRID.snow.Snow_i(GRID.snow.cT_domain_ub)=1./3.*GRID.snow.Snow_i(GRID.snow.cT_domain_ub+1);
            GRID.snow.Snow_w(GRID.snow.cT_domain_ub)=1./3.*GRID.snow.Snow_w(GRID.snow.cT_domain_ub+1);
            GRID.snow.Snow_a(GRID.snow.cT_domain_ub)=1./3.*GRID.snow.Snow_a(GRID.snow.cT_domain_ub+1);
            GRID.snow.Snow_i(GRID.snow.cT_domain_ub+1)=GRID.snow.Snow_i(GRID.snow.cT_domain_ub+1) - GRID.snow.Snow_i(GRID.snow.cT_domain_ub);
            GRID.snow.Snow_w(GRID.snow.cT_domain_ub+1)=GRID.snow.Snow_w(GRID.snow.cT_domain_ub+1) - GRID.snow.Snow_w(GRID.snow.cT_domain_ub);
            GRID.snow.Snow_a(GRID.snow.cT_domain_ub+1)=GRID.snow.Snow_a(GRID.snow.cT_domain_ub+1) - GRID.snow.Snow_a(GRID.snow.cT_domain_ub);
            T(GRID.snow.cT_domain_ub)=T(GRID.snow.cT_domain_ub+1);
            check_change=true;
        end

        if min(GRID.snow.Snow_i(GRID.snow.cT_domain))<=0.5.*PARA.technical.SWEperCell %avoid looping when unnecessary 
            for i=GRID.snow.cT_domain_ub:GRID.snow.cT_domain_lb-1  %check all snow cells except for the lowermost one for too small ice and water contents - merge with lower cell

                 if GRID.snow.Snow_i(i)<=0.5.*PARA.technical.SWEperCell

                     %------ modify snow and air grid -----------------------------
                     GRID.snow.cT_domain(GRID.snow.cT_domain_ub)=0;
                     GRID.snow.K_domain(GRID.snow.K_domain_ub)=0;
                     [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
                     [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

                     GRID.air.cT_domain(GRID.air.cT_domain_lb+1)=1;
                     GRID.air.K_domain(GRID.air.K_domain_lb+1)=1;
                     [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
                     [GRID.air.K_domain_lb, GRID.air.K_domain_ub] =   LayerIndex(GRID.air.K_domain);

                     %-------- rearrange SWE and T grids --------------------------        
                     GRID.snow.Snow_i(i+1)=GRID.snow.Snow_i(i+1)+GRID.snow.Snow_i(i);
                     GRID.snow.Snow_w(i+1)=GRID.snow.Snow_w(i+1)+GRID.snow.Snow_w(i);
                     GRID.snow.Snow_a(i+1)=GRID.snow.Snow_a(i+1)+GRID.snow.Snow_a(i);
                     GRID.snow.Snow_i(2:i)=GRID.snow.Snow_i(1:i-1);
                     GRID.snow.Snow_w(2:i)=GRID.snow.Snow_w(1:i-1);
                     GRID.snow.Snow_a(2:i)=GRID.snow.Snow_a(1:i-1);
                     T(i+1)=(T(i+1)+T(i))/2;
                     T(2:i)=T(1:i-1);
                 end
            end
            check_change=true;
        end

        if GRID.snow.Snow_i(GRID.snow.cT_domain_lb)<=0.5.*PARA.technical.SWEperCell && sum(GRID.snow.cT_domain)>=2 %lowermost grid cell has too little snow, but there still is 2 or more snow cells - merge with upper cell

            %------ modify snow and air grid ----------------------------------
            GRID.snow.cT_domain(GRID.snow.cT_domain_ub)=0;
            GRID.snow.K_domain(GRID.snow.cT_domain_ub)=0;
            [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
            [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

            GRID.air.cT_domain(GRID.air.cT_domain_lb+1)=1;
            GRID.air.K_domain(GRID.air.K_domain_lb+1)=1;
            [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
            [GRID.air.K_domain_lb, GRID.air.K_domain_ub] =   LayerIndex(GRID.air.K_domain);

            %-------- rearrange SWE and T grids --------------------------               
            GRID.snow.Snow_i(GRID.snow.cT_domain_lb)=GRID.snow.Snow_i(GRID.snow.cT_domain_lb)+GRID.snow.Snow_i(GRID.snow.cT_domain_lb-1);
            GRID.snow.Snow_w(GRID.snow.cT_domain_lb)=GRID.snow.Snow_w(GRID.snow.cT_domain_lb)+GRID.snow.Snow_w(GRID.snow.cT_domain_lb-1);
            GRID.snow.Snow_a(GRID.snow.cT_domain_lb)=GRID.snow.Snow_a(GRID.snow.cT_domain_lb)+GRID.snow.Snow_a(GRID.snow.cT_domain_lb-1);
            GRID.snow.Snow_i(2:GRID.snow.cT_domain_lb-1)=GRID.snow.Snow_i(1:GRID.snow.cT_domain_lb-2);
            GRID.snow.Snow_w(2:GRID.snow.cT_domain_lb-1)=GRID.snow.Snow_w(1:GRID.snow.cT_domain_lb-2);
            GRID.snow.Snow_a(2:GRID.snow.cT_domain_lb-1)=GRID.snow.Snow_a(1:GRID.snow.cT_domain_lb-2);

            T(GRID.snow.cT_domain_lb)=(T(GRID.snow.cT_domain_lb)+T(GRID.snow.cT_domain_lb-1))/2;
            T(2:GRID.snow.cT_domain_lb-1)=T(1:GRID.snow.cT_domain_lb-2);

            check_change=true;
        end


        if check_change==true %update grid spacings  
            
            assert( sum( sum( isnan(GRID.snow.Snow_i(GRID.snow.cT_domain) + GRID.snow.Snow_w(GRID.snow.cT_domain) + GRID.snow.Snow_a(GRID.snow.cT_domain)) ) ) == 0 , 'updateGRID_snow - error in Snow_i/a/w grid' );
            
            % snow grid
            soilTop = GRID.general.K_grid(GRID.snow.cT_domain_lb+1);  %tsvd in case of a lake, soilTop corresponds to lakeTop
            GRID.general.K_grid(GRID.snow.cT_domain)= soilTop - flipud(cumsum(flipud(GRID.snow.Snow_i(GRID.snow.cT_domain) + GRID.snow.Snow_w(GRID.snow.cT_domain) + GRID.snow.Snow_a(GRID.snow.cT_domain))));
            GRID.general.K_grid(GRID.air.cT_domain)=[GRID.general.K_grid(GRID.air.cT_domain_lb)+(-snowCellSize)*(GRID.air.cT_domain_lb-1):snowCellSize:GRID.general.K_grid(GRID.air.cT_domain_lb)]';
 
            assert( sum( isnan( GRID.general.K_grid ) )==0, 'updateGRID_snow - error in K grid after snow grid update')

            % air grid  zzz this is not in FLAKE version... jjj
            snowTop = GRID.general.K_grid(GRID.snow.cT_domain_ub);
            numAirCells = sum( GRID.air.cT_domain );
            GRID.general.K_grid(GRID.air.cT_domain) = flip( snowTop-snowCellSize:-snowCellSize:snowTop-snowCellSize*numAirCells );
            
            assert( sum( isnan( GRID.general.K_grid ) )==0, 'updateGRID_snow - error in K grid after air grid update')
            
            %GRID.general.K_grid(GRID.air.cT_domain) = [GRID.general.K_grid(GRID.air.cT_domain_lb)+(-2*snowCellSize)*(GRID.air.cT_domain_lb-1):2*snowCellSize:GRID.general.K_grid(GRID.air.cT_domain_lb)]';
        end 

        if (GRID.snow.Snow_i(GRID.snow.cT_domain_ub)<=0.5.*PARA.technical.SWEperCell && sum(GRID.snow.cT_domain)<2)  %remove last grid cell if snow threshold is reached

           % for water balance: add snow of last grid cell to runoff
           BALANCE.water.dr_snowmelt = BALANCE.water.dr_snowmelt - ( GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub) ).*1000;

           %------ modify snow and air grid ----------------------------------
           GRID.snow.cT_domain(GRID.snow.cT_domain_ub)=0;
           GRID.snow.K_domain(GRID.snow.cT_domain_ub)=0;
           [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
           [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);

           GRID.air.cT_domain(GRID.air.cT_domain_lb+1)=1;
           GRID.air.K_domain(GRID.air.K_domain_lb+1)=1;
           [GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
           [GRID.air.K_domain_lb, GRID.air.K_domain_ub] =   LayerIndex(GRID.air.K_domain);


           GRID.snow.Snow_i(GRID.air.cT_domain_lb)=0;
           GRID.snow.Snow_w(GRID.air.cT_domain_lb)=0;
           GRID.snow.Snow_a(GRID.air.cT_domain_lb)=0;
           T(GRID.air.cT_domain_lb)=0;

        end       


    end

    GRID.general.K_grid(GRID.air.cT_domain_lb)= GRID.general.K_grid(GRID.air.cT_domain_lb+1)-snowCellSize;
    GRID.general.K_grid(GRID.air.cT_domain_lb-1)= GRID.general.K_grid(GRID.air.cT_domain_lb+1)-2.*snowCellSize;

    GRID.general.cT_grid=( GRID.general.K_grid(1:end-1)+ GRID.general.K_grid(2:end))./2; %grid on which capacity and temperature information lives (midpoints of grid cells)
    GRID.general.cT_delta=(- GRID.general.cT_grid(1:end-1,1)+ GRID.general.cT_grid(2:end,1));
    GRID.general.K_delta=(- GRID.general.K_grid(1:end-1,1)+ GRID.general.K_grid(2:end,1));


    %bugfix as still situations may occur where K_delta<0  ... zzz also if lake level changes !?
    if ( sum( GRID.general.K_delta < 0 ) > 0 ) 
        disp('updateGRID_snow - bugfix K grid');
        %update grid spacings
        % snow grid
        if ~isempty(GRID.snow.cT_domain_ub) % snow cover
            GRID.general.K_grid(GRID.snow.cT_domain) = GRID.general.K_grid(GRID.snow.cT_domain_lb+1) - flipud(cumsum(flipud(GRID.snow.Snow_i(GRID.snow.cT_domain) + GRID.snow.Snow_w(GRID.snow.cT_domain) + GRID.snow.Snow_a(GRID.snow.cT_domain))));        
            surfaceTop = GRID.general.K_grid(GRID.snow.cT_domain_ub);
       else %lll no snow cover
%tsvd  LAKE cases added       surfaceTop = GRID.general.K_grid(GRID.soil.cT_domain_ub);
            if(~isempty(GRID.lake.ice.cT_domain_ub) ) % lake ice
                surfaceTop = GRID.general.K_grid(GRID.ice.water.cT_domain_ub);
            elseif(~isempty(GRID.lake.water.cT_domain_ub) ) % open lake
                surfaceTop = GRID.general.K_grid(GRID.lake.water.cT_domain_ub);
            else
                surfaceTop = GRID.general.K_grid(GRID.soil.cT_domain_ub);    % soil
            end
        end
        % air grid
        numAirCells = sum( GRID.air.cT_domain );
        GRID.general.K_grid(GRID.air.cT_domain) = flip( surfaceTop-snowCellSize:-snowCellSize:surfaceTop-snowCellSize*numAirCells );
        GRID.general.cT_grid  = (  GRID.general.K_grid(1:end-1)    + GRID.general.K_grid(2:end)    ) ./ 2;
        GRID.general.cT_delta = ( -GRID.general.cT_grid(1:end-1,1) + GRID.general.cT_grid(2:end,1) );
        GRID.general.K_delta  = ( -GRID.general.K_grid(1:end-1,1)  + GRID.general.K_grid(2:end,1)  );
    end
    
    assert( sum( isnan(T(GRID.snow.cT_domain)))==0, 'updateGRID_snow - error in T after grid update' );
    assert( sum( isnan(GRID.general.K_grid))==0, 'update_GRID_snow - error in Kgrid after grid update');
end
