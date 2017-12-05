function [ind_lb, ind_ub] = LayerIndex(vec) %#codegen
%clear all
%vec=[0 0 0 0 0 0 0 0 0 0 0 0 0 0]';


A=[1:length(vec)]';
C=A(logical(vec));

if isempty(C)==1
    ind_ub=[];
    ind_lb=[];
else
    ind_ub=C(1,1);
    ind_lb=C(end,1);
end