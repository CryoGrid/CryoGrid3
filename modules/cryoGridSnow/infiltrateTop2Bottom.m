function [snow_w, snow_a, waterFlux] = infiltrateTop2Bottom(snow_i, snow_w, snow_a, poreSpace, maxLiqWater, waterFlux)



for j=1:size(snow_w,1)
    delta_SWE=double((waterFlux - maxLiqWater(j))>0).*maxLiqWater(j) + double((waterFlux - maxLiqWater(j))<=0).*waterFlux;
    
    snow_w(j) = snow_w(j) + delta_SWE;
    if delta_SWE>=0
        
        snow_a(j) = snow_a(j) - delta_SWE;
        
    else
        snow_a(j)=(poreSpace(j).*snow_i(j) + poreSpace(j).*snow_w(j) - snow_w(j))./(1-poreSpace(j));
        
    end
    waterFlux=double((waterFlux - maxLiqWater(j))>0).*(waterFlux - maxLiqWater(j));
end

assert(sum(snow_a<0)==0,'Infiltrate Bottom2Top : Negative value of snow_a')
   