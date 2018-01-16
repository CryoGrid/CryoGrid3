function [GRID, PARA] = excessGroundIce(T, GRID, PARA)

if ~isempty(PARA.soil.mobileWaterDomain) && (sum(double(T(GRID.soil.cT_domain)>0 & GRID.soil.excessGroundIce==1))~=0) && isempty(GRID.snow.cT_domain_ub)
    disp('excess ice thawing');
    GRID.soil.excessGroundIce = GRID.soil.excessGroundIce==1 & T(GRID.soil.cT_domain)<=0;  %remove the thawed cell from the list
    [GRID meltwaterGroundIce PARA] = excessGroundIceThaw4(T, GRID, PARA);   %meltwaterGroundIce could be read out, but is not yet implemented
    GRID = updateGRID_excessice(GRID);
end