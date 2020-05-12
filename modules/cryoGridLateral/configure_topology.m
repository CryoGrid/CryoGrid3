function [ PARA ] = configure_topology( PARA, TOPO, SETUP )

% check if "TOPO" is of class "MesoTopology"
assert( isa( TOPO, 'MesoTopology'), 'TOPO must be an instance of the class MesoTopology' );

% here we create an "effective" micro-topology object which incorporates
% the topology on the micro- and meso-scales
effTopo = EffectiveMicroTopology(TOPO);

% introduce new variable to store the sizes of the meso/micro units
PARA.ensemble.microUnitSizes = TOPO.childSizes;

% finally, the topology/geometry information stored in the topology object
% are transferred to the "old" structure using PARA.ensemble

% the following choice of topological parameters is just an example for three connected linear tiles
PARA.ensemble.area = effTopo.areas;
if iscolumn(PARA.ensemble.area)
    PARA.ensemble.area=PARA.ensemble.area.';
end
PARA.ensemble.weight = PARA.ensemble.area;
PARA.ensemble.distanceBetweenPoints = effTopo.distances;

% adjacency of all lateral exchange networks
PARA.ensemble.adjacency_meso=effTopo.adjacency_meso;
PARA.ensemble.adjacency_micro=effTopo.adjacency_micro;

%water (surf/subsurf)
if SETUP.xWmeso
    PARA.ensemble.adjacency_water_surface=effTopo.adjacency;
    PARA.ensemble.adjacency_water_subsurface=effTopo.adjacency;
    PARA.ensemble.adjacency_water = max( PARA.ensemble.adjacency_water_surface, PARA.ensemble.adjacency_water_subsurface );
else
    PARA.ensemble.adjacency_water_surface=effTopo.adjacency_micro;
    PARA.ensemble.adjacency_water_subsurface=effTopo.adjacency_micro;
    PARA.ensemble.adjacency_water = max( PARA.ensemble.adjacency_water_surface, PARA.ensemble.adjacency_water_subsurface );
end

%heat
if SETUP.xHmeso
    PARA.ensemble.adjacency_heat=effTopo.adjacency;
else
    PARA.ensemble.adjacency_heat=effTopo.adjacency_micro;
end

%sediment
if SETUP.xEmeso
    PARA.ensemble.adjacency_sediment=effTopo.adjacency;
else
    PARA.ensemble.adjacency_sediment=effTopo.adjacency_micro;
end

%snow
PARA.ensemble.adjacency_snow=effTopo.adjacency_micro; % note that snow is distributed among all connected components of a micro unit

% parameters related to HEAT exchange
PARA.ensemble.thermal_contact_length = effTopo.contactLength.*PARA.ensemble.adjacency_heat;
PARA.ensemble.thermalDistance = PARA.ensemble.distanceBetweenPoints.*PARA.ensemble.adjacency_heat;

% parameters related to WATER exchange
PARA.ensemble.hydraulic_contact_length = effTopo.contactLength.*PARA.ensemble.adjacency_water;
PARA.ensemble.hydraulicDistance = PARA.ensemble.distanceBetweenPoints.*PARA.ensemble.adjacency_water;

% topographical relations
PARA.ensemble.initial_altitude = effTopo.altitudes; %in m a.s.l.
if iscolumn(PARA.ensemble.initial_altitude)
    PARA.ensemble.initial_altitude=PARA.ensemble.initial_altitude.';
end
PARA.ensemble.altitude = PARA.ensemble.initial_altitude;  
PARA.ensemble.surface_altitude = PARA.ensemble.initial_altitude;
PARA.ensemble.soil_altitude = PARA.ensemble.initial_altitude;

end