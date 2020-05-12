classdef EffectiveMicroTopology < Topology
    properties
        adjacency_micro(:,:) {mustBeInteger}
        adjacency_meso(:,:) {mustBeInteger}
        reservoir_altitude(1,1) {mustBeNumeric}
    end
    
    methods
        % constructor
        function obj = EffectiveMicroTopology(mesoTopo)
            if nargin==0
                super_args={};
            elseif nargin==1
                if isa(mesoTopo,'MesoTopology')
                    super_args{1}=sum(mesoTopo.childSizes);            
                else
                    error('Argument must be of type MesoTopology.');
                end
            else
                error('Too many input arguments.');
            end

            %ini
            obj@Topology(super_args{:});

            %post ini
            obj.type = 'effective';
            obj.scale = 'micro';
            
            % configure effective topology
            obj.setEffectiveAdjacency(mesoTopo);    % this function also sets the adjacency_meso/micro matrices
            obj.setEffectiveAltitudes(mesoTopo);
            obj.setEffectiveAreas(mesoTopo);
            obj.setEffectiveDistances(mesoTopo);
            obj.setEffectiveContactLength(mesoTopo);
            
            obj.reservoir_altitude=mesoTopo.reservoir_altitude;
            
        end    
        
        function obj = setEffectiveAdjacency(obj,parent)
            adjMicro=[];
            for i=1:parent.size
                adjMicro= blkdiag( adjMicro, parent.childUnits(i).adjacency);
            end
            adjMeso=zeros(sum(parent.childSizes));
            for i=1:parent.size
                for j=1:parent.size
                    if parent.adjacencyMicroIndices(i,j)
                        adjMeso( getGlobalIndex( i, parent.adjacencyMicroIndices(i,j), parent.childSizes ), ...
                                 getGlobalIndex( j, parent.adjacencyMicroIndices(j,i), parent.childSizes ) ) = 1;
                    end
                end
            end
            adjacencyFull = adjMicro + adjMeso;
            obj.adjacency_meso = adjMeso;
            obj.adjacency_micro = adjMicro;
            assert( issymmetric(adjacencyFull), 'adjacencyFull is not symmetric');
            obj.adjacency = adjacencyFull;
        end
        
        function obj = setEffectiveAltitudes(obj,parent)
            altitudesFull = zeros( parent.size, 1 );
            for i=1:parent.size
                for j=1:parent.childSizes(i)
                    altitudesFull( getGlobalIndex( i, j, parent.childSizes ) ) = parent.childUnits(i).altitudes(j);
                end
            end    
            obj.altitudes = altitudesFull;
        end
        
        function obj = setEffectiveAreas(obj,parent)
            areasFull = zeros( parent.size, 1 );
            for i=1:parent.size
                for j=1:parent.childSizes(i)
                    areasFull( getGlobalIndex( i, j, parent.childSizes ) ) = parent.childUnits(i).areas(j);
                end
            end    
            obj.areas = areasFull;
        end
        
        
        function obj = setEffectiveDistances(obj,parent)
            distMicro=[];
            for i=1:parent.size
                distMicro= blkdiag( distMicro, parent.childUnits(i).distances);
            end
            distMeso=zeros(sum(parent.childSizes));
            for i=1:parent.size
                for j=1:parent.size
                    if parent.adjacencyMicroIndices(i,j)
                        distMeso( getGlobalIndex( i, parent.adjacencyMicroIndices(i,j), parent.childSizes ), ...
                                 getGlobalIndex( j, parent.adjacencyMicroIndices(j,i), parent.childSizes ) ) = parent.distances(i,j);
                    end
                end
            end
            distFull = distMicro + distMeso;
            assert( issymmetric(distFull), 'distFull is not symmetric');
            obj.distances = distFull;
        end

        function obj = setEffectiveContactLength(obj,parent)
            contMicro=[];
            for i=1:parent.size
                contMicro= blkdiag( contMicro, parent.childUnits(i).contactLength);
            end
            contMeso=zeros(sum(parent.childSizes));
            for i=1:parent.size
                for j=1:parent.size
                    if parent.adjacencyMicroIndices(i,j)
                        contMeso( getGlobalIndex( i, parent.adjacencyMicroIndices(i,j), parent.childSizes ), ...
                                 getGlobalIndex( j, parent.adjacencyMicroIndices(j,i), parent.childSizes ) ) = ...
                                    parent.contactLength(i,j) .* parent.childUnits(i).area ./ parent.areas(i);      % here the correct scaling of micro-scale fluxes is ensured
                    end
                end
            end
            contFull = contMicro + contMeso;
            obj.contactLength = contFull;
        end       
        
        
        
    end
            
        
end