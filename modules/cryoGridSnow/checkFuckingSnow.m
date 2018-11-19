function [ GRID ] = checkFuckingSnow( GRID )
% Function to unbug the NaN bug in update grid snow because Snow_i or
% Snow_a get to have an NaN value

do_it=0;

if length(GRID.snow.Snow_i)-sum(GRID.snow.Snow_i==0)==1 && sum(isnan(GRID.snow.Snow_i))==1;
    GRID.snow.Snow_i(isnan(GRID.snow.Snow_i))=0;
    do_it=do_it+1;
end

if length(GRID.snow.Snow_a)-sum(GRID.snow.Snow_a==0)==1 && sum(isnan(GRID.snow.Snow_a))==1;
    GRID.snow.Snow_a(isnan(GRID.snow.Snow_a))=0;
    do_it=do_it+1;
end

if do_it==2;
    GRID.snow.cT_domain = logical(0.*GRID.snow.cT_domain);
    GRID.snow.K_domain  = logical(0.*GRID.snow.K_domain) ;
    
    [GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
    [GRID.snow.K_domain_lb, GRID.snow.K_domain_ub] =   LayerIndex(GRID.snow.K_domain);
end

end

