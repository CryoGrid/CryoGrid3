function p = satPresWater(T)
%Magnus Formula
p=0.622.* 6.112 .* 100 .* exp(17.62.*(T-273.15)./(243.12-273.15+T));