function [loop_PARA loopI_PERMUT] = generateLOOPVAR(LOOPVAR, LI)


%example
% LI=5;
% LOOPVAR.ice_albedo = [0.1 0.2 0.3];
% LOOPVAR.ice_extinction = [7.0 8.4 9.0];
% LOOPVAR.water_extinction = [0.2 0.5 0.8];
% LOOPVAR.FLAKE_fetch = [50 150 250];

var_names = fieldnames(LOOPVAR);

for i=1:length(var_names)    
    var_names{i,2}=['i_' num2str(i)];
    var_names{i,3}=num2str(length(eval(['LOOPVAR.' var_names{1}])));   
end

start_str=[];
mind_str=[];
end_str=[];
for i=1:length(var_names)
    start_str=[start_str 'for ' var_names{i,2} '=1:' var_names{i,3} '; '];
    mind_str=[mind_str var_names{i,2} ' '];
    end_str=[end_str 'end; '];   
end
loop_index=[];
eval([start_str 'loop_index=[loop_index; ' mind_str ']; ' end_str])

for i=1:length(var_names)
    eval(['loopI_PERMUT.' var_names{i,1} '=' 'loop_index(:,' num2str(i) ');'])
    eval(['loop_PARA.' var_names{i,1} '=' 'LOOPVAR.' var_names{i,1} '(loop_index(' num2str(LI) ',' num2str(i) '));'  ])
    
    
end











