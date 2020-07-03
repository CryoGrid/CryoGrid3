function [precond_water] = checkPreconditionWaterExchange( T, GRID )
  
    precond_water = 0; %tsvd IS
	precondition_waterExchange = double( T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub) ); % matches with conditions of infiltration
    
% 	if(labindex<numlabs)
%         labSend( precondition_waterExchange, labindex+1, 22);
%     end
% %tsvd IS
%     labBarrier(); % ddd ccc
%     if(T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub))  % only update precond_water if water exchange for current tile is activated (i.e. precondition_waterExchange is true)
%         if(labindex>1)
%             precondition_waterExchange = precondition_waterExchange + labReceive( labindex-1, 22);
%         end
% %tsvd IS        precond_water = precondition_waterExchange > 1;	% at least two workers need matching conditions
%         precond_water = precondition_waterExchange == 2;	% at least two workers need matching conditions
%     end  
    
%     j_prevWorker = mod(labindex - 2, numlabs) + 1;  j_nextWorker = mod(labindex, numlabs) + 1;
%     precondition_waterExchange_j = labSendReceive(j_nextWorker,j_prevWorker,precondition_waterExchange,22)
%     
%     precond_water = precondition_waterExchange + precondition_waterExchange_j > 1;	% at least two workers need matching conditions
    
    j_prevWorker = mod(labindex - 2, numlabs) + 1; 
    j_nextWorker = mod(labindex, numlabs) + 1; 
    precondition_waterExchange_prev = labSendReceive(j_nextWorker,j_prevWorker,precondition_waterExchange,22); 
    precondition_waterExchange_next = labSendReceive(j_prevWorker,j_nextWorker,precondition_waterExchange,22);
    
    precond_water = precondition_waterExchange + precondition_waterExchange_prev +  precondition_waterExchange_next > 1;	% at least two workers need matching conditions
    

%%%%%   Original:  %%%%%%%%%

% function [precond_water] = checkPreconditionWaterExchange( T, GRID )
% 
% 	precondition_waterExchange = double( T(GRID.soil.cT_domain_ub)>0 && isempty(GRID.snow.cT_domain_ub) ); % matches with conditions of infiltration
% 
% 	for j=1:numlabs
%         if j~=labindex
%             labSend( precondition_waterExchange, j, 2);
%         end
%     end
% 	
%     for j=1:numlabs
%         if j~=labindex
%             precondition_waterExchange = precondition_waterExchange + labReceive( j, 2);
%         end
%     end
%             
% 	precond_water = precondition_waterExchange > 1;	% at least two workers need matching conditions
    