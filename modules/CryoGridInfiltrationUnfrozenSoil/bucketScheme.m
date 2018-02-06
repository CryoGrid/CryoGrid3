function [wc, surface_runoff, lacking_water]=bucketScheme(T, wc, dwc_dt, GRID, PARA, external_flux)

T=T(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
porosity=1-GRID.soil.cT_mineral-GRID.soil.cT_organic; 

% A=sum(K_delta(1:30).*wc(1:30)+dwc_dt(1:30))
lacking_water=0;
i=1;
i_max=70;  % maximum infiltration depth, must be defined somehow before
while  T(i)>0 && i<=i_max
    max_water=K_delta(i).*PARA.soil.fieldCapacity;  %maximum amount of water (in m) that a grid cell can hold
    min_water=K_delta(i).*PARA.soil.residualWC;     %minimum amount of water which stays in a cell
    
    actual_water= max( min_water,  wc(i).*K_delta(i)+dwc_dt(i) );       % should be dwc (already multiplied with timestep) %JAN: this violates the WB
    
    lacking_water = lacking_water + (wc(i).*K_delta(i)+dwc_dt(i)<min_water) .* (wc(i).*K_delta(i)+dwc_dt(i)-min_water);
    
    dwc_dt(i+1)=dwc_dt(i+1) + max(0, actual_water-max_water);  %when excess water, move it to next grid cell
    wc(i)=min(max_water, actual_water)./K_delta(i); 
       
    i=i+1;
end

excess_water=dwc_dt(i)+external_flux; %add external flux

% B=sum(K_delta(1:30).*wc(1:30))+excess_water
i=i-1;

while i>=1 && excess_water>0 
    max_water=K_delta(i).*porosity(i);
    actual_water=wc(i).*K_delta(i)+excess_water;
    wc(i)=min(actual_water, max_water)./K_delta(i);
    excess_water=max(0, actual_water-wc(i).*K_delta(i));
       
    i=i-1;
end


surface_runoff=excess_water;

% C=sum(K_delta(1:30).*wc(1:30))+surface_runoff

