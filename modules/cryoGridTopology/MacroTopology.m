classdef MacroTopology < Topology

    properties
        childUnits(:,1) MesoTopology
    end
    
    methods
        function obj = MacroTopology(arg1,arg2)
            if nargin==0
                super_args={};
            elseif nargin==1
                super_args{1}=arg1;
            elseif nargin>1
                super_args{1}=arg1;
                super_args{2}=arg2;
            elseif nargin>2
                    error('Specify valid number of input arguments.');
            end
            
            % ini
            obj@Topology(super_args{:});
            
            obj.scale='macro';
        end
    end
    
end

