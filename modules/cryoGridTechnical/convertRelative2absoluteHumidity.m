function q=convertRelative2absoluteHumidity(FORCING)
p=1005*100;

q=(double((FORCING.data.Tair)<=0).*(FORCING.data.RH.*satPresIce(FORCING.data.Tair+273.15)) + double(FORCING.data.Tair>0).*(FORCING.data.RH.*satPresWater(FORCING.data.Tair+273.15)))./p;