function [ geomSetup ] = iLoadGeomSetup( filename, config )
% Functon to load the geomSetup

geomSetup=load(filename);
geomSetup=geomSetup.geomSetup(config);

end

