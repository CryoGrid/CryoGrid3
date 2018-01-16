function [GRID, PARA, wc, meltwaterGroundIce] = excessGroundIceInfiltration(T, wc, GRID, PARA)
meltwaterGroundIce=0;
if ~isempty(PARA.soil.mobileWaterDomain) && (sum(double(T(GRID.soil.cT_domain)>0 & GRID.soil.excessGroundIce==1))~=0) && isempty(GRID.snow.cT_domain_ub)
    disp('excess ice thawing');
    GRID.soil.excessGroundIce = GRID.soil.excessGroundIce==1 & T(GRID.soil.cT_domain)<=0;   %remove the thawed cell from the list
    [GRID, meltwaterGroundIce, wc] = excessGroundIceThaw4Infiltration(T, wc, GRID, PARA);   %meltwaterGroundIce could be read out, but is not yet implemented
    
    %[GRID, wc] = updateGRID_excessiceInfiltration(wc, GRID);
    %JAN: updateGRID not necessary as long as no distinct water domain, the
    %regridding of soil/air domain happens already in the Thaw4 function
    
    % modification due to infiltration --> finally change cT_water to wc
    %GRID.soil.cT_water = wc;
    % JAN: do NOT update cT_water as this stores the latest frozen water
    % content and is used to check whether the LUT needs to be updated in
    % the infiltration module
    
end