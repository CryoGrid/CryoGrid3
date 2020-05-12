classdef Topology < handle
    properties
        type = 'none'
        scale = 'none' 
        size(1,1) {mustBeInteger} = 1
        adjacency(:,:) {mustBeInteger} = 0
        altitudes(:,1) {mustBeNumeric} = 0
        areas(:,1) {mustBeNumeric} = 1
        area(1,1) {mustBeNumeric} = 1
        distances(:,:) {mustBeNumeric} = 0
        contactLength(:,:) {mustBeNumeric} = 0
    end

    methods
        % CONTRUCTORS
        %       Pre-initialization —    Compute arguments for superclass constructors.
        %       Object initialization — Call superclass constructors.
        %       Post initialization —   Perform any operations related to the subclass,
        %                               including referencing and assigning to the object, 
        %                               call class methods, passing the object to functions, and so on.
        
        function obj = Topology(arg1,arg2,distance)
            % pre ini
            ...
            % ini
            ...
            % post ini
            if nargin>0
                if nargin==1
                    if isnumeric(arg1)
                        size=arg1;
                        type='simple';
                    elseif ischar(arg1)
                        size=1;
                        type=arg1;
                    end
                elseif nargin>1
                    type=arg1;
                    if isscalar(arg2)
                        size=arg2;
                    else
                        size=1;
                    end
                end
                obj.type = type;
                obj.size = size;
                obj.adjacency = zeros(size);
                obj.altitudes = zeros(size,1);
                obj.areas = ones(size,1);
                obj.update_area();
                obj.distances = zeros(size);
                obj.contactLength = zeros(size);
                
                % check special types
                if strcmp(type,'translational')
                    obj.adjacency = diag( ones(size-1,1), 1 ) + diag( ones(size-1,1), -1 );
                    if nargin>2
                        L=1;
                        obj.areas = L.*distance.*ones(size,1);
                        obj.distances = distance .* obj.adjacency;
                        obj.contactLength = L .* obj.adjacency;
                        obj.update_area();
                    end
                elseif strcmp(arg1,'radial')
                    if nargin==2
                        error('Third argument [distance] must be specified for radial geometry.');
                    elseif nargin>2
                        obj.adjacency = diag( ones(size-1,1), 1 ) + diag( ones(size-1,1), -1 );
                        obj.areas = 2.*pi.*distance.^2 .* ( cumsum(ones(size,1))-1 ) ;
                        obj.areas(1) = pi ./ 4 .* distance.^2;
                        obj.distances = distance .* obj.adjacency;
                        obj.contactLength = distance .* pi .* ( diag( (2*cumsum(ones(size-1,1))-1), 1 ) + diag( (2*cumsum(ones(size-1,1))-1), -1 ) );
                        obj.update_area()
                    end
                end
                %error('If more than one argument is passed, the first must specify the type.\n Valid types: "radial", "translational", "simple"');

            else

                ...
            end
            
        end
        
        % GETTER and SETTER functions
        function set.areas(obj,vector)
            obj.validateVectorDimension(vector);
            obj.areas=vector;
            obj.update_area();
        end
        
        function set.altitudes(obj,vector)
            obj.validateVectorDimension(vector);
            obj.altitudes=vector;
        end
        
        function set.adjacency(obj,matrix)
            obj.validateMatrixDimension(matrix);
            obj.adjacency=matrix;
        end
        
        function set.distances(obj,matrix)
            obj.validateMatrixDimension(matrix);
            obj.distances=matrix;
        end
        
        function set.contactLength(obj,matrix)
            obj.validateMatrixDimension(matrix);
            obj.contactLength=matrix;
        end
       
        % class functions
        function update_area(obj)
            obj.area=sum(obj.areas);
        end
        
        function validateMatrixDimension(obj,matrix)
            dim=size(matrix);
            assert( dim(1)==obj.size, 'Matrix dimension does not match size of object.');
            assert( dim(2)==obj.size, 'Matrix dimension does not match size of object.');
        end
        
        function validateVectorDimension(obj,vector)
            assert( length(vector)==obj.size, 'Vector dimension does not match size of object.');
        end
            
        
    end
    
end

