function GRID = createStratigraphy(PARA, GRID)

    i=1;   
    soilParam=[-100 PARA.soil.layer_properties(i, 2:6) ];
    i=2;
    while i<=size(PARA.soil.layer_properties,1)
        soilParam=[soilParam; PARA.soil.layer_properties(i, 1)-0.01 PARA.soil.layer_properties(i-1, 2:6)];
        soilParam=[soilParam; PARA.soil.layer_properties(i, 1)+0.01 PARA.soil.layer_properties(i, 2:6)];
        i=i+1;
    end
    soilParam=[soilParam; 1e6 PARA.soil.layer_properties(i-1, 2:6)];


    %interpolate to grid
    GRID.soil.cT_water = interp1(soilParam(:,1),soilParam(:,2),GRID.general.cT_grid(GRID.soil.cT_domain),'linear');
    GRID.soil.cT_mineral = interp1(soilParam(:,1),soilParam(:,3),GRID.general.cT_grid(GRID.soil.cT_domain),'linear');
    GRID.soil.cT_organic = interp1(soilParam(:,1),soilParam(:,4),GRID.general.cT_grid(GRID.soil.cT_domain),'linear');
    GRID.soil.cT_soilType = interp1(soilParam(:,1),soilParam(:,5),GRID.general.cT_grid(GRID.soil.cT_domain),'nearest');
    GRID.soil.cT_natPor = interp1(soilParam(:,1),soilParam(:,6),GRID.general.cT_grid(GRID.soil.cT_domain),'linear');

    GRID.soil.K_water = interp1(soilParam(:,1),soilParam(:,2),GRID.general.K_grid(GRID.soil.K_domain),'linear');
    GRID.soil.K_mineral = interp1(soilParam(:,1),soilParam(:,3),GRID.general.K_grid(GRID.soil.K_domain),'linear');
    GRID.soil.K_organic = interp1(soilParam(:,1),soilParam(:,4),GRID.general.K_grid(GRID.soil.K_domain),'linear');
    GRID.soil.K_soilType = interp1(soilParam(:,1),soilParam(:,5),GRID.general.K_grid(GRID.soil.K_domain),'nearest');

end