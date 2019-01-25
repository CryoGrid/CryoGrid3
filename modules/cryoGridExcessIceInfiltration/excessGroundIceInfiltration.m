function [GRID, PARA, wc, meltwaterGroundIce] = excessGroundIceInfiltration(T, wc, GRID, PARA)

meltwaterGroundIce=0;

if ~isempty(PARA.soil.mobileWaterDomain) && (sum(double(T(GRID.soil.cT_domain)>0 & GRID.soil.excessGroundIce==1))~=0) && isempty(GRID.snow.cT_domain_ub)
   
    GRID.soil.excessGroundIce = GRID.soil.excessGroundIce==1 & T(GRID.soil.cT_domain)<=0;   %remove the thawed cell from the list
    
    % Now distinguish cases. It the soil is frozen or not between the
    % mobile water and the surface.
    status = excessGroundIceOverlyingConfig( T, wc, GRID );
    
    if strcmp(status,'unfrozen');
         disp('excessGroundIceInfiltration - excess ice thawing');
        [GRID, meltwaterGroundIce, wc] = excessGroundIceThaw4Infiltration(T, wc, GRID, PARA);   % meltwaterGroundIce could be read out, but is not yet implemented
        
    else
         disp('excessGroundIceInfiltration - excess ice thawing - no water ponding');
        [GRID, meltwaterGroundIce, wc] = excessGroundIceThawLateral(T, wc, GRID, PARA);
       
    end
    
end