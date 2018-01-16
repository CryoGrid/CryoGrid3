function [cT_water, cT_mineral, cT_organic, K_water, K_mineral, K_organic, K_grid]=...
   removeWater(T, cT_water, cT_mineral, cT_organic, cT_natPor, K_water, K_mineral, K_organic, K_delta, K_grid, targetSaturation)
   

%function[cT_water, cT_mineral, cT_organic,  K_grid]=...
 %   removeWater(T, cT_water, cT_mineral, cT_organic, cT_natPor, K_delta, K_grid, targetSaturation)

    superSaturatedCells=find(cT_water(:,1)>cT_natPor(:,1) & T(:,1)>0);
    for i=superSaturatedCells
        
        newGridCellSize=(cT_mineral(i,1)+cT_organic(i,1)).*K_delta(i,1)./(1-cT_natPor(i,1));   %in [m]
        K_grid(i+1:end)=K_grid(i+1:end)-(K_delta(i,1)-newGridCellSize);
        cT_water(i,1)=cT_natPor(i,1).*targetSaturation;
        cT_mineral(i,1)=cT_mineral(i,1).*K_delta(i,1)./newGridCellSize;
        cT_organic(i,1)=cT_organic(i,1).*K_delta(i,1)./newGridCellSize;
      % mobileWater=(cT_water(i,1)-cT_natPor(i,1)).*K_delta(i,1); %in [m]
       
      % K_grid(i+1:end)=K_grid(i+1:end)-mobileWater;
       
       
       
        K_water(i,1)=(cT_water(i-1,1)+cT_water(i,1))/2;
        K_mineral(i,1)=(cT_mineral(i-1,1)+cT_mineral(i,1))/2;
        K_organic(i,1)=(cT_organic(i-1,1)+cT_organic(i,1))/2;
            
        K_water(i+1,1)=(cT_water(i+1,1)+cT_water(i,1))/2;
        K_mineral(i+1,1)=(cT_mineral(i+1,1)+cT_mineral(i,1))/2;
        K_organic(i+1,1)=(cT_organic(i+1,1)+cT_organic(i,1))/2;
      
    
    end
    
    
    %mineral+org content (cT_mineral(i,1)+cT_organic(i,1)).*K_delta(i,1)
    
    %soll (1-cT_natPor).*newK_delta ergeben