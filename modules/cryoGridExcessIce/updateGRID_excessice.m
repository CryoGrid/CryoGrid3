function GRID = updateGRID_excessice(GRID)
%--- update soil and water grid after subsidence ----------------------

% change soil with 100% water to water cell
soilGRIDsize = sum(GRID.soil.cT_domain);

<<<<<<< HEAD

% JAN: set all cells above without soil to air (water runs off and no snow cover
% present)
% this needs to be improved!!!
=======
% set all cells above without soil to air (water runs off and no snow cover
% present)
>>>>>>> origin/xice_mpi_polygon_TC
GRID.air.cT_domain(GRID.soil.cT_domain) = (GRID.soil.cT_organic==0 & GRID.soil.cT_mineral==0);
GRID.air.K_domain(GRID.soil.K_domain) = (GRID.soil.K_organic==0 & GRID.soil.K_mineral==0);

GRID.soil.cT_domain(GRID.soil.cT_domain) = (GRID.soil.cT_organic>0 | GRID.soil.cT_mineral>0);
GRID.soil.K_domain(GRID.soil.K_domain)   = (GRID.soil.K_organic>0 | GRID.soil.K_mineral>0);

<<<<<<< HEAD


=======
>>>>>>> origin/xice_mpi_polygon_TC
if soilGRIDsize ~= sum(GRID.soil.cT_domain)
    
    disp('subsidence - updating grid information');
    
    %water_bucket = GRID.soil.cT_water(1);
    cT_no_water = (GRID.soil.cT_organic>0 | GRID.soil.cT_mineral>0);
    K_no_water  = (GRID.soil.K_organic>0 | GRID.soil.K_mineral>0);
    
    [GRID.soil.cT_domain_lb GRID.soil.cT_domain_ub] = LayerIndex(GRID.soil.cT_domain);
    [GRID.soil.K_domain_lb GRID.soil.K_domain_ub]   = LayerIndex(GRID.soil.K_domain);
    
    [GRID.air.cT_domain_lb GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
    [GRID.air.K_domain_lb GRID.air.K_domain_ub]   = LayerIndex(GRID.air.K_domain);
<<<<<<< HEAD

    
%     GRID.water.cT_domain(max([GRID.air.cT_domain_lb+1 GRID.ice.cT_domain_lb+1]) : GRID.soil.cT_domain_ub-1) = 1;
%     GRID.water.K_domain(max([GRID.air.K_domain_lb+1 GRID.ice.K_domain_lb+1]) : GRID.soil.K_domain_ub-1) = 1;
%     
%     [GRID.water.cT_domain_lb GRID.water.cT_domain_ub] = LayerIndex(GRID.water.cT_domain);
%     [GRID.water.K_domain_lb GRID.water.K_domain_ub]   = LayerIndex(GRID.water.K_domain);
    
    %-- update all other soil grid infos if size has changed
    


=======
    
    %-- update all other soil grid infos if size has changed
>>>>>>> origin/xice_mpi_polygon_TC
    % adjust cT grid fields
    GRID.soil.cT_water = GRID.soil.cT_water(cT_no_water);
    GRID.soil.cT_mineral = GRID.soil.cT_mineral(cT_no_water);
    GRID.soil.cT_organic = GRID.soil.cT_organic(cT_no_water);
    GRID.soil.cT_soilType = GRID.soil.cT_soilType(cT_no_water);
    GRID.soil.cT_natPor = GRID.soil.cT_natPor(cT_no_water);
    GRID.soil.excessGroundIce = GRID.soil.excessGroundIce(cT_no_water);
    GRID.soil.conductivity = GRID.soil.conductivity(cT_no_water, :);
    GRID.soil.capacity = GRID.soil.capacity(cT_no_water, :);
    GRID.soil.liquidWaterContent = GRID.soil.liquidWaterContent(cT_no_water, :);
    GRID.soil.cT_frozen = GRID.soil.cT_frozen(cT_no_water);
    GRID.soil.cT_thawed = GRID.soil.cT_thawed(cT_no_water);
    GRID.soil.K_frozen = GRID.soil.cT_frozen; % this is ok since K_frozen and cT_frozen
    GRID.soil.K_thawed = GRID.soil.cT_thawed; % are the same from the start (see initialize.m)
    
    % adjust K grid fields
    GRID.soil.soilGrid = GRID.soil.soilGrid(K_no_water);
    GRID.soil.K_water = GRID.soil.K_water(K_no_water);
    GRID.soil.K_mineral = GRID.soil.K_mineral(K_no_water);
    GRID.soil.K_organic = GRID.soil.K_organic(K_no_water);
    GRID.soil.K_soilType = GRID.soil.K_soilType(K_no_water);
    
<<<<<<< HEAD
    %     s = fieldnames(GRID.soil);
%     for i=1:length(s)
%         if isempty(strfind(char(s(i)),'domain')) && isempty(strfind(char(s(i)),'K_frozen')) && isempty(strfind(char(s(i)),'K_thawed'))% exclude all with name 'domain'
%             if isempty(strfind(char(s(i)),'K_')) && isempty(strfind(char(s(i)),'soilGrid'))
%                 %evaluate on cT grid
%                 eval(['GRID.soil.' char(s(i)) '=' 'GRID.soil.' char(s(i)) '(cT_no_water,:);']);
%             else
%                 %evaluate on K grid
%                 eval(['GRID.soil.' char(s(i)) '=' 'GRID.soil.' char(s(i)) '(K_no_water,:);']);
%             end
%         end
%     end
=======
>>>>>>> origin/xice_mpi_polygon_TC
end