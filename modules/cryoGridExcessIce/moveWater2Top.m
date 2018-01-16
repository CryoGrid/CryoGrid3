function [cT_water, cT_mineral, cT_organic, K_water, K_mineral, K_organic]=moveWater2Top(T, cT_water, cT_mineral, cT_organic, cT_natPor, K_water, K_mineral, K_organic, K_delta, cT_firstGroundCell)
    superSaturatedCells=find(cT_water(:,1)>cT_natPor(:,1) & T(:,1)>0);
    cellsChanged=superSaturatedCells;
    for i=1:size(superSaturatedCells,1)
       mobileWater=(cT_water(superSaturatedCells(i),1)-cT_natPor(superSaturatedCells(i),1)).*K_delta(superSaturatedCells(i),1); %in [m]
       j=cT_firstGroundCell;
       while mobileWater>0 && j<=size(cT_water,1)
          waterAdded=min(mobileWater, K_delta(j,1).*(1-cT_water(j,1)));
          cT_water(j,1)=cT_water(j,1)+waterAdded./K_delta(j,1);
          cT_water(superSaturatedCells(i),1)= cT_water(superSaturatedCells(i),1)-waterAdded./K_delta(superSaturatedCells(i),1);
          if cT_organic(j,1)+cT_mineral(j,1)>0
            organicSubtracted=waterAdded.*cT_organic(j,1)./(cT_organic(j,1)+cT_mineral(j,1));
            mineralSubtracted=waterAdded.*cT_mineral(j,1)./(cT_organic(j,1)+cT_mineral(j,1));
          else
            organicSubtracted=0;
            mineralSubtracted=0;
          end
          cT_organic(j,1)=cT_organic(j,1)-organicSubtracted./K_delta(j,1);
          cT_organic(superSaturatedCells(i),1)= cT_organic(superSaturatedCells(i),1)+organicSubtracted./K_delta(superSaturatedCells(i),1);
          
          cT_mineral(j,1)=cT_mineral(j,1)-mineralSubtracted./K_delta(j,1);
          cT_mineral(superSaturatedCells(i),1)= cT_mineral(superSaturatedCells(i),1)+mineralSubtracted./K_delta(superSaturatedCells(i),1);
          
          mobileWater=mobileWater-waterAdded;
          cellsChanged=[cellsChanged; j];
          j=j+1;
       end
    end
    
    for i=cellsChanged'
       if i==cT_firstGroundCell
           K_water(i,1)=cT_water(i,1);
           K_mineral(i,1)=cT_mineral(i,1);
           K_organic(i,1)=cT_organic(i,1);
       else
           K_water(i,1)=(cT_water(i-1,1)+cT_water(i,1))/2;
           K_mineral(i,1)=(cT_mineral(i-1,1)+cT_mineral(i,1))/2;
           K_organic(i,1)=(cT_organic(i-1,1)+cT_organic(i,1))/2;
           
           K_water(i+1,1)=(cT_water(i+1,1)+cT_water(i,1))/2;
           K_mineral(i+1,1)=(cT_mineral(i+1,1)+cT_mineral(i,1))/2;
           K_organic(i+1,1)=(cT_organic(i+1,1)+cT_organic(i,1))/2;
       end
    end