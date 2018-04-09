function [wc, surface_runoff, lacking_water]=bucketScheme(T, wc, dwc_dt, GRID, PARA, external_flux)

T=T(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
porosity=1-GRID.soil.cT_mineral-GRID.soil.cT_organic;
soilType = GRID.soil.cT_soilType;
[~,i_max]=min(abs((PARA.location.altitude - GRID.soil.soilGrid)-PARA.soil.infiltration_limit_altitude));

fieldCapacity = zeros(size(soilType));
residualWaterContent = zeros(size(soilType));
for i=1:size(PARA.soil.soilTypes,1)
	fieldCapacity(soilType==i) = PARA.soil.soilTypes( i, 2 );
	residualWaterContent(soilType==i) = PARA.soil.soilTypes( i, 1 );
end

lacking_water=0;
i=1;

while  T(i)>0 && i<=i_max
    max_water=K_delta(i).*fieldCapacity(i);  %maximum amount of water (in m) that a grid cell can hold
    min_water=K_delta(i).*residualWaterContent(i);     %minimum amount of water which stays in a cell (independent of soil type, but should be if "freezing = drying")
    
    actual_water= max( min_water,  wc(i).*K_delta(i)+dwc_dt(i) );       % should be dwc (already multiplied with timestep) %JAN: this violates the WB
    
    lacking_water = lacking_water + (actual_water - (wc(i).*K_delta(i)+dwc_dt(i) ) );

    dwc_dt(i+1)=dwc_dt(i+1) + max(0, actual_water-max_water);  %when excess water, move it to next grid cell
    wc(i)=min(max_water, actual_water)./K_delta(i);
    i=i+1;
end

excess_water=dwc_dt(i)+external_flux; %add external flux
excess_water=excess_water - lacking_water; % remove potential mismatches (e.g. when evaporation in cell with low water content)

lacking_water = -(excess_water<0)*excess_water;  % this accounts for violations of the water balance
i=i-1;

while i>=1 && excess_water>0
    max_water=K_delta(i).*porosity(i);
    actual_water=wc(i).*K_delta(i)+excess_water;
    wc(i)=min(actual_water, max_water)./K_delta(i);
    excess_water=max(0, actual_water-wc(i).*K_delta(i));
    i=i-1;
end

surface_runoff=(excess_water>0)*excess_water; % surface runoff only if excess_water>0