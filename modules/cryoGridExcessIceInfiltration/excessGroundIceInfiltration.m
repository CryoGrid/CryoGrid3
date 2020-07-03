function [GRID, PARA, wc, meltwaterGroundIce] = excessGroundIceInfiltration(T, wc, GRID, PARA)
meltwaterGroundIce=0;
if ~isempty(PARA.soil.mobileWaterDomain) && (sum(double(T(GRID.soil.cT_domain)>0 & GRID.soil.excessGroundIce==1))~=0) && isempty(GRID.snow.cT_domain_ub)
    disp('check 2 index')
    labindex
    disp('excessGroundIceInfiltration - excess ice thawing');
    GRID.soil.excessGroundIce = GRID.soil.excessGroundIce==1 & T(GRID.soil.cT_domain)<=0;   %remove the thawed cell from the list
    [GRID, meltwaterGroundIce, wc] = excessGroundIceThaw4Infiltration(T, wc, GRID, PARA);   %meltwaterGroundIce could be read out, but is not yet implemented
        
end