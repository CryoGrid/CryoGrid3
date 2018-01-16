function k_old=conductivityFreeWater(k_old, T, soilWater)

k_freeWater=5; %set conductivity for mobile water to 5
i=1;

while T(i)>0 && soilWater(i)>=1
    k_old(i)=k_freeWater;
    i=i+1;
end