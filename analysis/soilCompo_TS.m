function [ TS ] = soilCompo_TS( CellArray, compo_name, ind_depth, disp )
% create a time serie of one component of the soil from output

if strcmp(compo_name,'water')
    ind_compo=1;
elseif strcmp(compo_name,'mineral')
    ind_compo=2;
elseif strcmp(compo_name,'organic')
    ind_compo=3;
else
    error('soil component "%s" does not exit',compo_name)
end

myper=@(x) permute(x,[2 3 1]); % Create an anonymous function handle with the appropriate permutation
B=cellfun(myper,CellArray,'UniformOutput',false); % use it with cellfun
C=cell2mat(B); % convert into mat
D=permute(C,[3 2 1]); % Permute so that lines = depths, rows = time, layers=water / minerale / organic

TS=D(ind_depth,:,ind_compo)'; % extract the good TS

if disp==1;
    plot(TS)
end

end