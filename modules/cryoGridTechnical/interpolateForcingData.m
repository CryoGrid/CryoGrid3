function [FORCING]= interpolateForcingData(t, FORCING)

assert(imag(t)==0,'interpolateForcingData Error : t is complex')

posit=floor((t-FORCING.data.t_span(1,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1)))+1;
  
FORCING.i.snowfall=FORCING.data.snowfall(posit,1)+(FORCING.data.snowfall(posit+1,1)-FORCING.data.snowfall(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.rainfall=FORCING.data.rainfall(posit,1)+(FORCING.data.rainfall(posit+1,1)-FORCING.data.rainfall(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.Lin=FORCING.data.Lin(posit,1)+(FORCING.data.Lin(posit+1,1)-FORCING.data.Lin(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.Sin=FORCING.data.Sin(posit,1)+(FORCING.data.Sin(posit+1,1)-FORCING.data.Sin(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.Tair=FORCING.data.Tair(posit,1)+(FORCING.data.Tair(posit+1,1)-FORCING.data.Tair(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.wind=FORCING.data.wind(posit,1)+(FORCING.data.wind(posit+1,1)-FORCING.data.wind(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.q=FORCING.data.q(posit,1)+(FORCING.data.q(posit+1,1)-FORCING.data.q(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));

FORCING.i.p=FORCING.data.p(posit,1)+(FORCING.data.p(posit+1,1)-FORCING.data.p(posit,1)).*(t-FORCING.data.t_span(posit,1))./(FORCING.data.t_span(2,1)-FORCING.data.t_span(1,1));