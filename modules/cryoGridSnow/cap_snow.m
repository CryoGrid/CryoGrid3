function cap_snow=cap_snow(Snow_i, Snow_w, Snow_a, PARA)

c_i = PARA.constants.c_i; %1.9*10^6;%[J/mK]
c_w = PARA.constants.c_w; %4.2*10^6; %[J/mK]

cap_snow = (Snow_i.* c_i + Snow_w.* c_w)./ (Snow_i + Snow_w + Snow_a);
cap_snow(find(isnan(cap_snow(:,1))==1),1)=0.3.*c_i;

