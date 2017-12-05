function cond_snow = cond_snow(snow_i, snow_w, snow_a)

total=snow_w + snow_i + snow_a;


%cond_snow=conductivity2(snow_w./total, snow_i./total, 0, 0, 2);
cond_snow = 2.2.*((snow_w+snow_i)./total).^1.88; % Yen (1981)

cond_snow(find(isnan(cond_snow)==1),1)=0.3;
    