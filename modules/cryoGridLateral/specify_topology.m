function [ TOPO ] = specify_topology( SETUP )
      
    TOPO=MesoTopology(SETUP.mesoType,SETUP.Nmeso,SETUP.Dmeso,SETUP.altitude,SETUP.slope);
    
    for i=1:SETUP.Nmeso
        TOPO.childUnits(i)=MicroTopology(TOPO,i,SETUP.microType,SETUP.microArea,SETUP.microElev);
    end
    
    % specify adjacency relations on the micro-level
    if strcmp(SETUP.microType,'polygonCRT')
        TOPO.adjacencyMicroIndices = 3 .* TOPO.adjacency;
    else
        TOPO.adjacencyMicroIndices = TOPO.adjacency;    %connect micro-tiles with index 1
    end
    
end