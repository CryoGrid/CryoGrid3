function altitude = getAltitude( PARA, GRID )
   %zzz same as get soil altitude!  used?
%tsvd account for lake surface (snow can accumulate on lake ice)
if ~isempty( GRID.lake.water.cT_domain_ub )% lake water is present
    altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.lake.water.cT_domain_ub);
elseif ~isempty( GRID.lake.ice.cT_domain_ub ) % lake ice is present
    altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.lake.ice.cT_domain_ub); 
    
else
    altitude = PARA.location.initial_altitude - GRID.general.K_grid(GRID.soil.cT_domain_ub); % no lake
end

end