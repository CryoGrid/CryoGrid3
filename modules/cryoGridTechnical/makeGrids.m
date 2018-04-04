function GRID=makeGrids(PARA)

GRID.snow.snowCellSize=PARA.technical.SWEperCell/(PARA.snow.rho_snow/PARA.constants.rho_w);
GRID.snow.snowGrid=[-1.*(PARA.technical.maxSWE./(PARA.technical.SWEperCell)+2).*GRID.snow.snowCellSize:GRID.snow.snowCellSize:-GRID.snow.snowCellSize]';
%tsvd added 
GRID.lake.water.waterGrid=[0:0.02:PARA.water.depth-0.02]';  %with water (it is recommended to use a grid resolution of about 2cm)
%GRID.lake.water.waterGrid=[]'; %no water
if ~isempty(GRID.lake.water.waterGrid)
    GRID.soil.soilGrid=PARA.technical.subsurfaceGrid + (2*GRID.lake.water.waterGrid(end)-GRID.lake.water.waterGrid(end-1));
else
GRID.soil.soilGrid=PARA.technical.subsurfaceGrid;
end

%tsvd   K_grid =[GRID.snow.snowGrid;  GRID.soil.soilGrid]; %grid on which the conductivty information lives (edges of grid cells)
K_grid =[GRID.snow.snowGrid; GRID.lake.water.waterGrid; GRID.soil.soilGrid]; %grid on which the conductivty information lives (edges of grid cells)
cT_grid=(K_grid(1:end-1)+K_grid(2:end))/2; %grid on which heat capacity and temperature information lives (midpoints of grid cells)
cT_delta=(-cT_grid(1:end-1,1)+cT_grid(2:end,1));
K_delta=(-K_grid(1:end-1,1)+K_grid(2:end,1));

GRID.general.cT_grid=cT_grid;
GRID.general.K_grid=K_grid;
GRID.general.cT_delta=cT_delta;
GRID.general.K_delta=K_delta;

%set air grid
GRID.air.cT_domain= false(size(GRID.general.cT_grid));
GRID.air.cT_domain(1:length(GRID.snow.snowGrid))=1;
GRID.air.K_domain= false(size(GRID.general.K_grid));
GRID.air.K_domain(1:length(GRID.snow.snowGrid))=1;
[GRID.air.cT_domain_lb, GRID.air.cT_domain_ub] = LayerIndex(GRID.air.cT_domain);
[GRID.air.K_domain_lb, GRID.air.K_domain_ub]   = LayerIndex(GRID.air.K_domain);

%set soil grid
GRID.soil.cT_domain= false(size(GRID.general.cT_grid));
%tsvd    GRID.soil.cT_domain(GRID.general.cT_grid>0)=1;
GRID.soil.cT_domain(end-length(GRID.soil.soilGrid)+2:end)=1; %required if there is a lake on top 
[GRID.soil.cT_domain_lb, GRID.soil.cT_domain_ub] = LayerIndex(GRID.soil.cT_domain);
GRID.soil.K_domain= false(size(GRID.general.K_grid));
GRID.soil.K_domain(GRID.soil.cT_domain_ub:end)=1;
[GRID.soil.K_domain_lb, GRID.soil.K_domain_ub]   = LayerIndex(GRID.soil.K_domain);

%set snow grid
GRID.snow.cT_domain=false(size(GRID.general.cT_grid));
GRID.snow.K_domain=false(size(GRID.general.K_grid));
[GRID.snow.cT_domain_lb, GRID.snow.cT_domain_ub] = LayerIndex(GRID.snow.cT_domain);
[GRID.snow.K_domain_lb, GRID.snow.K_domain_ub]   = LayerIndex(GRID.snow.K_domain);

% JAN: currently these are "secondary" domains which belong to the
% "primary" soil domain. Lake cells are treated like soil cells during
% integration of the heat conduction, but the surface properties are
% different and mixing occurs in summer. Initially, the lake domain is
% empty, but it is updated in the first time step.

%tsvd replaced by Flake conditions below   - zzz Jan, if you update your model version with using GRID.lake.water.* instead of GRID.lake.*, our handling of GRID.lake is consistent 
% %set water and ice cover grid
% GRID.lake.cT_domain = false(size(GRID.general.cT_grid));
% %GRID.lake.cT_domain(GRID.air.cT_domain_lb+1:GRID.soil.cT_domain_ub-1)=1;
% GRID.lake.K_domain= false(size(GRID.general.K_grid));
% %GRID.lake.K_domain(GRID.air.K_domain_lb+1:GRID.soil.K_domain_ub-1)=1;
% [GRID.lake.cT_domain_lb, GRID.lake.cT_domain_ub] = LayerIndex(GRID.lake.cT_domain);
% [GRID.lake.K_domain_lb, GRID.lake.K_domain_ub]   = LayerIndex(GRID.lake.K_domain);
% 
% GRID.lake.water.cT_domain= false(size(GRID.general.cT_grid));
% GRID.lake.water.K_domain= false(size(GRID.general.K_grid));
% [GRID.lake.water.cT_domain_lb, GRID.lake.water.cT_domain_ub] = LayerIndex(GRID.lake.water.cT_domain);
% [GRID.lake.water.K_domain_lb, GRID.lake.water.K_domain_ub]   = LayerIndex(GRID.lake.water.K_domain);
% 
% GRID.lake.ice.cT_domain=false(size(GRID.general.cT_grid));
% GRID.lake.ice.K_domain=false(size(GRID.general.K_grid));                
% [GRID.lake.ice.cT_domain_lb, GRID.lake.ice.cT_domain_ub] = LayerIndex(GRID.lake.ice.cT_domain);
% [GRID.lake.ice.K_domain_lb, GRID.lake.ice.K_domain_ub]   = LayerIndex(GRID.lake.ice.K_domain); 

%tsvd 
%set water and ice cover grid
GRID.lake.water.cT_domain= false(size(GRID.general.cT_grid));
GRID.lake.water.cT_domain(GRID.air.cT_domain_lb+1:GRID.soil.cT_domain_ub-1)=1;
GRID.lake.water.K_domain= false(size(GRID.general.cT_grid));
GRID.lake.water.K_domain(GRID.air.K_domain_lb+1:GRID.soil.K_domain_ub-1)=1;
[GRID.lake.water.cT_domain_lb GRID.lake.water.cT_domain_ub] = LayerIndex(GRID.lake.water.cT_domain);
[GRID.lake.water.K_domain_lb GRID.lake.water.K_domain_ub]   = LayerIndex(GRID.lake.water.K_domain);
GRID.lake.ice.cT_domain=logical(GRID.air.cT_domain.*0);
GRID.lake.ice.K_domain=logical(GRID.air.K_domain.*0);                  
[GRID.lake.ice.cT_domain_lb GRID.lake.ice.cT_domain_ub] = LayerIndex(GRID.lake.ice.cT_domain);
[GRID.lake.ice.K_domain_lb GRID.lake.ice.K_domain_ub]   = LayerIndex(GRID.lake.ice.K_domain); 
GRID.lake.ice.z_ice = 0;
GRID.lake.ice.dz_dt_ice = 0;
GRID.lake.ice.dE_dt_cond_residual=0;
GRID.lake.ice.melt_flag=false;