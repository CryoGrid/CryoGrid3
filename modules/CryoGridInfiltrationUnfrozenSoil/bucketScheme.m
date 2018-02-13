function [wc, surface_runoff, lacking_water]=bucketScheme(T, wc, dwc_dt, GRID, PARA, external_flux,index)

T=T(GRID.soil.cT_domain);
K_delta=GRID.general.K_delta(GRID.soil.cT_domain);  %in m
porosity=1-GRID.soil.cT_mineral-GRID.soil.cT_organic;
soilType = GRID.soil.cT_soilType;
if numlabs < 2;
    [~,i_max]=min(abs((PARA.location.altitude(index)-GRID.soil.soilGrid)-PARA.soil.alt_infiltration_limit));
else
    [~,i_max]=min(abs((PARA.ensemble.altitude(index)-GRID.soil.soilGrid)-PARA.ensemble.alt_infiltration_limit));
end

% to be changed!
fieldCapacity = zeros(size(soilType));
residualWaterContent = zeros(size(soilType));
for i=1:size(PARA.soil.soilTypes,1)
    fieldCapacity(soilType==i) = PARA.soil.soilTypes( i, 2 );
    residualWaterContent(soilType==1) = PARA.soil.soilTypes( i, 1 );
end

lacking_water=0;
i=1;
% i_max=length(dwc_dt)-1; % 200;  % maximum infiltration depth, must be defined somehow before, includes also water body on to of soil

while  T(i)>0 && i<=i_max
    
    max_water=K_delta(i).*fieldCapacity(i);  %maximum amount of water (in m) that a grid cell can hold
    min_water=K_delta(i).*residualWaterContent(i);     %minimum amount of water which stays in a cell (independent of soil type, but should be if "freezing = drying")
    
    actual_water= max( min_water,  wc(i).*K_delta(i)+dwc_dt(i) );       % should be dwc (already multiplied with timestep) %JAN: this violates the WB
    
    lacking_water = lacking_water + (wc(i).*K_delta(i)+dwc_dt(i)-actual_water);
    
    dwc_dt(i+1)=dwc_dt(i+1) + max(0, actual_water-max_water);  %when excess water, move it to next grid cell
    wc(i)=min(max_water, actual_water)./K_delta(i);
    
    i=i+1;
end

excess_water=dwc_dt(i)+external_flux; %add external flux
lacking_water = lacking_water + (excess_water<0)*excess_water;  % this accounts for violations of the water balance

i=i-1;

% if index==1;
%     fprintf('Inside')
%     wc
%     excess_water
%     external_flux
%     K_delta
% end

while i>=1 && excess_water>0
    max_water=K_delta(i).*porosity(i);
    actual_water=wc(i).*K_delta(i)+excess_water;
    wc(i)=min(actual_water, max_water)./K_delta(i);
    excess_water=max(0, actual_water-wc(i).*K_delta(i));
    
    i=i-1;
end


surface_runoff=(excess_water>0)*excess_water; % surface runoff only if excess_water>0

