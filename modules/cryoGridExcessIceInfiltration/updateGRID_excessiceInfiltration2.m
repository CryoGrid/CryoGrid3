function [GRID] = updateGRID_excessiceInfiltration2(meltwaterGroundIce, GRID)

    % pass excess meltwater to storage variable
    GRID.lake.residualWater = meltwaterGroundIce;

    % update GRID domains of water body
    if GRID.soil.cT_organic(1)+GRID.soil.cT_mineral(1)<=1e-6    % upper soil cell pure air/water
        % general water body extent
        cT_waterBody = GRID.soil.cT_organic+GRID.soil.cT_mineral<=1e-6;
        GRID.lake.cT_domain(logical(GRID.air.cT_domain+GRID.snow.cT_domain)) = 0;
        GRID.lake.cT_domain(GRID.soil.cT_domain) = cT_waterBody;
        [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
         GRID.lake.K_domain(logical(GRID.air.K_domain+GRID.snow.K_domain)) = 0;
         GRID.lake.K_domain(GRID.lake.cT_domain_ub:GRID.lake.cT_domain_lb+1) = 1;
         GRID.lake.K_domain(GRID.lake.cT_domain_lb+2:end) = 0;
         [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub] = LayerIndex(GRID.lake.K_domain);

%         % distinction frozen and unfrozen parts
%         unfrozenWaterBody = GRID.soil.cT_organic+GRID.soil.cT_mineral<=1e-6 & T(GRID.soil.cT_domain)>0;
%         frozenWaterBody = GRID.soil.cT_organic+GRID.soil.cT_mineral<=1e-6 & T(GRID.soil.cT_domain)<=0;
%         GRID.lake.water.cT_domain = GRID.lake.cT_domain & T>0;
%         [GRID.lake.water.cT_domain_lb, GRID.lake.water.cT_domain_ub] = LayerIndex(GRID.lake.water.cT_domain);
%         GRID.lake.ice.cT_domain = GRID.lake.cT_domain & T<=0;
%         [GRID.lake.ice.cT_domain_lb, GRID.lake.ice.cT_domain_ub] = LayerIndex(GRID.lake.ice.cT_domain); %these might be two domains
%         % K domains not implemented so far
          
    else
        GRID.lake.cT_domain = false(size(GRID.general.cT_grid));
        GRID.lake.K_domain = false(size(GRID.general.K_grid));
        [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
        [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub] = LayerIndex(GRID.lake.K_domain);
    end

end