classdef MicroTopology < Topology
    
    properties
        parentUnit(:,1) MesoTopology
        parentIndex(1,1) {mustBeInteger}
    end
    
    methods
        % CONTRUCTORS
        %       Pre-initialization —    Compute arguments for superclass constructors.
        %       Object initialization — Call superclass constructors.
        %       Post initialization —   Perform any operations related to the subclass,
        %                               including referencing and assigning to the object, 
        %                               call class methods, passing the object to functions, and so on.
        
        function obj = MicroTopology(parent,parentIndex,arg1,arg2,arg3)
            

            if nargin<3
                super_args={};
            elseif nargin==3
                super_args{1}=arg1;
            elseif nargin>3
                super_args{1}=arg1;
                super_args{2}=arg2;
                if nargin==5
                    super_args{3}=arg3; % distance between tiles
                end
                if nargin>5
                    error('Specify valid number of input arguments (1-5).');
                end
            end

            % ini
            obj@Topology(super_args{:});

            % post-ini
            obj.scale='micro';

            if nargin>1
                % link to parent object and back
                obj.parentUnit = parent;
                obj.parentIndex = parentIndex;
                assert( parentIndex <= parent.size, 'Parent index out of range.');
                obj.parentUnit.childUnits(obj.parentIndex)=obj;
            end

            % special cases
            if nargin>2
                if strcmp(arg1,'polygonCRT')
                    if nargin==3
                        areasCRT=140.*[0.3,0.6,0.1];    % default values [Nitzbon et al. (2019)]
                        elevationCRT=[0.0,0.4,0.3];     % default values [Nitzbon et al. (2019)]
                    else
                        areasCRT = arg2;
                        elevationCRT = arg3;
                    end

                    obj.type = arg1;
                    obj.size = 3;

                    obj.altitudes = obj.parentUnit.altitudes(obj.parentIndex) + elevationCRT;     %add altitude of containing MESO-struct.
                    obj.areas = areasCRT;

                    R_C=sqrt(obj.areas(1)./pi);
                    R_R=sqrt((obj.areas(1)+obj.areas(2))./pi);
                    R_T=sqrt(sum(obj.areas)./pi);
                    D_CR=(R_C+R_R)./2;
                    D_RT=R_T-R_R./2-R_C./2;
                    L_CR=2*pi*R_C;
                    L_RT=2*pi*R_R;
                    obj.distances = [[0,D_CR,0];[D_CR,0,D_RT];[0,D_RT,0]];
                    obj.contactLength = [[0,L_CR,0];[L_CR,0,L_RT];[0,L_RT,0]];
                    obj.adjacency = obj.distances > 0;
                    obj.update_area();
                elseif strcmp(arg1,'simple')
                    obj.type = arg1;
                    obj.size = 1;
                    obj.altitudes = obj.parentUnit.altitudes(obj.parentIndex);     %add altitude of containing MESO-struct.
                end                    
            end
        end
    end  
end

