function [wc, GRID] = updateGRID_erosion( PARA, GRID, wc )



while (GRID.soil.cT_mineral(1)+GRID.soil.cT_organic(1)+wc(1)<1e-6)
    disp('infiltration - update GRID - removing air cell')

    % adjust air and soil domains and boundaries
    GRID.air.cT_domain(GRID.soil.cT_domain_ub)=1;
    GRID.air.K_domain(GRID.soil.K_domain_ub)=1;
    GRID.air.cT_domain_lb=GRID.air.cT_domain_lb+1;
    GRID.air.K_domain_lb=GRID.air.K_domain_lb+1;
    GRID.soil.cT_domain(GRID.soil.cT_domain_ub)=0;
    GRID.soil.K_domain(GRID.soil.K_domain_ub)=0;
    GRID.soil.cT_domain_ub=GRID.soil.cT_domain_ub+1;
    GRID.soil.K_domain_ub=GRID.soil.K_domain_ub+1;
    GRID.soil.soilGrid(1)=[];

    wc(1)=[];

    GRID.soil.cT_organic(1)=[];
    GRID.soil.cT_natPor(1)=[];
    GRID.soil.cT_actPor(1)=[];
    GRID.soil.cT_mineral(1)=[];
    GRID.soil.cT_soilType(1)=[];

    GRID.soil.excessGroundIce(1)=[];

end







% update GRID water body / lake domain
if GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-9    % upper soil cell pure air/water
    cT_waterBody = GRID.soil.cT_organic+GRID.soil.cT_mineral<=1e-9;
    GRID.lake.cT_domain(logical(GRID.air.cT_domain+GRID.snow.cT_domain)) = 0;
    GRID.lake.cT_domain(GRID.soil.cT_domain) = cT_waterBody;
    [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
    GRID.lake.K_domain(logical(GRID.air.K_domain+GRID.snow.K_domain)) = 0;
    GRID.lake.K_domain(GRID.lake.cT_domain_ub:GRID.lake.cT_domain_lb+1) = 1;
    GRID.lake.K_domain(GRID.lake.cT_domain_lb+2:end) = 0;
    [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub] = LayerIndex(GRID.lake.K_domain);
else
    GRID.lake.cT_domain = false(size(GRID.general.cT_grid));
    GRID.lake.K_domain = false(size(GRID.general.K_grid));
    [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
    [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub] = LayerIndex(GRID.lake.K_domain);
end

% update LUT (to be sure, in any case when fluxes where applied to
% the current realization)
fprintf('\t\t\t\t\t reinitializing LUT...');
GRID = initializeSoilThermalProperties(GRID, PARA);