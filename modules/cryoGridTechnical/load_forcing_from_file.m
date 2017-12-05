function [FORCING, success] = load_forcing_from_file(PARA)

load(['forcing/' PARA.forcing.filename])

FORCING.data.rainfall=FORCING.data.rainfall.*PARA.forcing.rain_fraction;
FORCING.data.snowfall=FORCING.data.snowfall.*PARA.forcing.snow_fraction;

if std(FORCING.data.t_span(2:end,1)-FORCING.data.t_span(1:end-1,1))>1e-8
    disp('timestamp of forcing data is not in regular intervals -> check, fix and restart')
    success=0;
    return
end

%here, consistency checks, RH->q calculation, set threhsolds for wind, etc. could be placed

FORCING.data.wind(FORCING.data.wind<0.5)=0.5; %set min wind speed to 0.5 m/sec to avoid breakdown of turbulence



%set pressure to mean pressure at corresponding altitude (international
%altitude formula) if now provided by the forcing time series
if ~isfield(FORCING.data, 'p')
    FORCING.data.p=FORCING.data.Tair.*0 + 1013.25.*100.*(1-0.0065./288.15.*PARA.location.altitude).^5.255;
end

success=1;

%initialize
FORCING.i.snowfall=0;
FORCING.i.rainfall=0;
FORCING.i.Lin=0;
FORCING.i.Sin=0;
FORCING.i.Tair=0;
FORCING.i.wind=0;
FORCING.i.RH=0;
FORCING.i.q=0;
FORCING.i.p=0;