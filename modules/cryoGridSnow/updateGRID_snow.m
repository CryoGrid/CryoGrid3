function [GRID, T] = updateGRID_snow(T, GRID, PARA)

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
        GRID.snow.Snow_a(GRID.snow.cT_domain_ub) = (GRID.snow.SWEinitial./(PARA.snow.rho_snow./1000) - GRID.snow.SWEinitial);      
        GRID.snow.SWEinitial=0;       
        
        % -------- update K grid -------------------------------------
        GRID.general.K_grid(GRID.snow.K_domain_ub)=-1.*( GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + GRID.snow.Snow_a(GRID.snow.cT_domain_ub) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub) );
        T(GRID.snow.cT_domain_ub)=T(GRID.air.cT_domain_lb);
       
    end
    
else   %snow exists

    check_change=false;      
    GRID.general.K_grid(GRID.snow.cT_domain_ub) = GRID.general.K_grid(GRID.snow.cT_domain_ub+1) -...
        ( GRID.snow.Snow_i(GRID.snow.cT_domain_ub) + GRID.snow.Snow_w(GRID.snow.cT_domain_ub) + GRID.snow.Snow_a(GRID.snow.cT_domain_ub)); %updates the position of the uppermost snow grid cell

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
      GRID.general.K_grid(GRID.snow.cT_domain)=GRID.general.K_grid(GRID.snow.cT_domain_lb+1) - flipud(cumsum(flipud(GRID.snow.Snow_i(GRID.snow.cT_domain) + GRID.snow.Snow_w(GRID.snow.cT_domain) + GRID.snow.Snow_a(GRID.snow.cT_domain))));
      GRID.general.K_grid(GRID.air.cT_domain)=[GRID.general.K_grid(GRID.air.cT_domain_lb)+(-snowCellSize)*(GRID.air.cT_domain_lb-1):snowCellSize:GRID.general.K_grid(GRID.air.cT_domain_lb)]';
    end 
        
    if (GRID.snow.Snow_i(GRID.snow.cT_domain_ub)<=0.5.*PARA.technical.SWEperCell && sum(GRID.snow.cT_domain)<2)  %remove last grid cell if snow threshold is reached
          
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

% update grid spacings
GRID.general.K_grid(GRID.air.cT_domain_lb)= GRID.general.K_grid(GRID.air.cT_domain_lb+1)-snowCellSize;
GRID.general.K_grid(GRID.air.cT_domain_lb-1)= GRID.general.K_grid(GRID.air.cT_domain_lb+1)-2.*snowCellSize;
GRID.general.cT_grid=( GRID.general.K_grid(1:end-1)+ GRID.general.K_grid(2:end))./2; %grid on which capacity and temperature information lives (midpoints of grid cells)
GRID.general.cT_delta=(- GRID.general.cT_grid(1:end-1,1)+ GRID.general.cT_grid(2:end,1));
GRID.general.K_delta=(- GRID.general.K_grid(1:end-1,1)+ GRID.general.K_grid(2:end,1));