function [snow_w, snow_a, waterFlux] = infiltrateBottom2Top(snow_i, snow_w, snow_a, waterFlux)
  
j=size(snow_w, 1);


while waterFlux>0 && j>=1
    delta_SWE = double((waterFlux - snow_a(j))>0).*snow_a(j) + double((waterFlux - snow_a(j))<=0).*waterFlux;
    snow_w(j) = snow_w(j) + delta_SWE;
    snow_a(j) = snow_a(j) - delta_SWE;
    waterFlux=waterFlux - delta_SWE;
    j=j-1;
end

assert(sum(snow_a<0)==0,'Infiltrate Bottom2Top : Negative value of snow_a')


