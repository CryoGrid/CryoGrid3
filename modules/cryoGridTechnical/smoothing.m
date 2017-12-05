function [Lin_n, Sin_n, Tair_n, wind_n, RH_n, rainfall_n, snowfall_n]=smoothing(Lin, Sin, Tair, wind, RH, rainfall, snowfall, t_span)

step=t_span(2,1)-t_span(1,1);

smoothing=(24/6)^-1;

step=round(smoothing/step/2);
Lin_n=Lin;
Sin_n=Sin;
Tair_n=Tair;
wind_n=wind;
RH_n=RH;
rainfall_n=rainfall;
snowfall_n=snowfall;


for i=step:size(t_span,1)-step-1
    Lin_n(i,1)=nanmean(Lin(i-step+1:i+step,1));
    Sin_n(i,1)=nanmean(Sin(i-step+1:i+step,1));
    Tair_n(i,1)=nanmean(Tair(i-step+1:i+step,1));
    wind_n(i,1)=nanmean(wind_n(i-step+1:i+step,1));
    RH_n(i,1)=nanmean(RH_n(i-step+1:i+step,1));
    snowfall_n(i,1)=nanmean(snowfall(i-step+1:i+step,1));
    rainfall_n(i,1)=nanmean(rainfall(i-step+1:i+step,1));
end
    