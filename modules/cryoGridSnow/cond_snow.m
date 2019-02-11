function cond_snow = cond_snow(snow_i, snow_w, snow_a)

total=snow_w + snow_i + snow_a;

assert(sum(imag(snow_w))==0,'cond_snow : snow_w is complex')
assert(sum(imag(snow_i))==0,'cond_snow : snow_i is complex')
assert(sum(imag(snow_a))==0,'cond_snow : snow_a is complex')


%cond_snow=conductivity2(snow_w./total, snow_i./total, 0, 0, 2);
cond_snow = 2.2.*((snow_w+snow_i)./total).^1.88; % Yen (1981)

assert(sum(imag(cond_snow))==0,'cond_snow : cond_snow is complex at flag 1')

cond_snow(find(isnan(cond_snow)==1),1)=0.3;

assert(sum(imag(cond_snow))==0,'cond_snow : cond_snow is complex at flag 2')

    