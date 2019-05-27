function [ PARA, waterflux ] = calculateLateralWaterOscillations( PARA, waterflux )
% Check if the lateral fluxes generate oscillations in the water content.
% So identify alternating patterns and brakes it.

pattern=PARA.ensemble.lastWaterChange;
pattern=circshift(pattern,[0,-1]);
pattern(end)=waterflux;

if sum(isnan(pattern))==0
    pos=sum(find(pattern>0));
    if mod(pos,2)==0 && pos>0 && pos<6
        waterflux=waterflux/2;
        pattern(end)=waterflux;
        fprintf('\t\t\t calculateLateralWaterOscillations : dividing water flux\n')
    end
end

PARA.ensemble.lastWaterChange=pattern;

end