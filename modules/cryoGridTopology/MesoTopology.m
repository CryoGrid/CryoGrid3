classdef MesoTopology < Topology

    properties
        parentUnit(1,1) MacroTopology
        parentIndex(1,1) {mustBeInteger}
        childUnits(:,1) MicroTopology
        adjacencyMicroIndices(:,:) {mustBeInteger}
        reservoir_altitude(1,1) {mustBeNumeric}
    end
    properties (SetAccess=private)
        childSizes(:,1) {mustBeInteger}
    end
    
    methods
        function obj = MesoTopology(arg1,arg2,distance,altitude,slope)
            if nargin==0
                super_args={};
            elseif nargin==1
                super_args{1}=arg1;
            elseif nargin>1
                super_args{1}=arg1;
                super_args{2}=arg2;
                if nargin>2
                    super_args{3}=distance;
                end
            elseif nargin>5
                error('Specify valid number of input arguments.');
            end
            
            % ini
            obj@Topology(super_args{:});
            
            % post-ini
            obj.scale='meso';
            obj.reservoir_altitude = 0;
            
            % initialize childUnits and link to this meso object
            A(1,obj.size)=MicroTopology();
            obj.childUnits = A;
            for i=1:obj.size
                obj.childUnits(i).parentUnit = obj;
                obj.childUnits(i).parentIndex = i;
            end
            
            % link to parent % so far only implemeted for one meso unit
            obj.parentIndex=1;
            obj.parentUnit.childUnits(obj.parentIndex)=obj;
            
            % special cases
            % SLOPE
            if nargin==4
                obj.altitudes = obj.altitudes + altitude;
            elseif nargin==5
                obj.altitudes = obj.altitudes + altitude + ( cumsum(ones(obj.size,1))-1 ) .* distance .* slope;
            end         
            
        end
        
        % SETTER
        function set.childUnits(obj,child)
            obj.childUnits=child;
            obj.update_childSizes();
        end
        
        function set.adjacencyMicroIndices(obj,matrix)
            obj.validateMatrixDimension(matrix);
            for i=1:obj.size
                assert( max( matrix(i,:) ) <= obj.childUnits(i).size, 'Adjacency indices do not match size of micro units.');
            end
            obj.adjacencyMicroIndices=matrix;
        end
        
        % UPDATER
        function update_childSizes(obj)
           obj.childSizes = [ obj.childUnits.size ];            
        end        
    end
end

