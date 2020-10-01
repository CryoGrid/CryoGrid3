function [zL1Tpos_ub,zL1Tpos_lb,zL2Tpos_ub,zL2Tpos_lb,zL3Tpos_ub,zL3Tpos_lb] = find_thawedLayers(zsoil,Z,diagTL)
% find depths of thawed layers or talik 

% Z: input matrix (Tsoil or LWC/WC)
% zsoil : layer mid-points
% diagTL determines how diagnostic is calculated:  0: calculation of thawed layers based on T>0 criterium, 1: talik occurence is calc based on LWC>threshold criterium

fid=fopen('test.txt','wt');

clear ind* Dind* zL*
SwitchPlot=0;
%zL1Tpos_ub(size(Z,2))=NaN; zL1Tpos_lb(size(Z,2))=NaN; zL2Tpos_ub(size(Z,2))=NaN; zL2Tpos_lb(size(Z,2))=NaN; zL3Tpos_ub(size(Z,2))=NaN; zL3Tpos_lb(size(Z,2))=NaN;
zL1Tpos_ub=nan(1,length(Z));zL1Tpos_lb=nan(1,length(Z));  zL2Tpos_ub=nan(1,length(Z));zL2Tpos_lb=nan(1,length(Z));  zL3Tpos_ub=nan(1,length(Z));zL3Tpos_lb=nan(1,length(Z));  

for ti=1:size(Z,2) % time loop
%for ti=1:2 % time loop
    if(diagTL==0) % diagnose thaw layers based on T>0
        indThaw=find(Z(:,ti)>0);
    elseif(diagTL==1) % diagnose thaw layers based on LWC>threshold
        indThaw=find(sum(Z,2)>0.9*size(Z,2)); % Z = LWC./WC  - if condition fullfilled, than full thaw (LWC/WC=1) for all time steps per year
        %indThaw=find(Z(:,ti)>1.0); % talik based on LWC.... (more than 50% of pore water must be liquid)
    end
    DindThaw=diff(indThaw);
    DindThaw_gt1=find(DindThaw>1);

%    if isempty(indThaw) % no thaw

%     elseif length(indThaw)==1 % only 1 layer of positive temperatures
%         zL1Tpos_ub(ti)=zsoil(indThaw);
%         zL1Tpos_lb(ti)=zsoil(indThaw);
    if ~isempty(indThaw) % at least one layer has thawing conditions
        zL1Tpos_ub(ti)=zsoil(indThaw(1));
        switch length(DindThaw_gt1) 
            case 0 % only 1 thaw layer
%                fprintf(fid,' 1 thaw layer')
                zL1Tpos_lb(ti)=zsoil(indThaw(end));
                zL2Tpos_ub(ti)=NaN; zL2Tpos_lb(ti)=NaN;
                zL3Tpos_ub(ti)=NaN; zL3Tpos_lb(ti)=NaN;
            case 1 % 2 thaw layers
%                fprintf(fid,' 2 thaw layers')
                zL1Tpos_lb(ti)=zsoil(indThaw(DindThaw_gt1));
                zL2Tpos_ub(ti)=zsoil(indThaw(DindThaw_gt1+1));
                zL2Tpos_lb(ti)=zsoil(indThaw(end));   
                zL3Tpos_ub(ti)=NaN; zL3Tpos_lb(ti)=NaN;
            case 2 % 3 thaw layers
%                fprintf(fid,' 3 thaw layers')
                zL1Tpos_lb(ti)=zsoil(indThaw(DindThaw_gt1(1)));
                zL2Tpos_ub(ti)=zsoil(indThaw(DindThaw_gt1(1)+1));
                zL2Tpos_lb(ti)=zsoil(indThaw(DindThaw_gt1(2)));
                zL3Tpos_ub(ti)=zsoil(indThaw(DindThaw_gt1(2)+1));
                zL3Tpos_lb(ti)=zsoil(indThaw(end));
            otherwise
%                fprintf(fid,' >3 thaw layers')
               if(length(DindThaw_gt1>2)); disp(['WARNING: number of taliks is: ',num2str(length(DindThaw_gt1))]); end
               zL1Tpos_lb(ti)=zsoil(indThaw(DindThaw_gt1(1)));
               zL2Tpos_ub(ti)=zsoil(indThaw(DindThaw_gt1(1)+1));
               zL2Tpos_lb(ti)=zsoil(indThaw(DindThaw_gt1(2)));
               zL3Tpos_ub(ti)=zsoil(indThaw(DindThaw_gt1(2)+1));
               zL3Tpos_lb(ti)=zsoil(indThaw(3));
        end
    end    
end
%zLminTpos=min([zL1Tpos_lb,zL2Tpos_lb,zL3Tpos_lb]);

if(SwitchPlot)
    time=datetime(ts,'ConvertFrom','datenum');

    figure(1) % plot Thaw Layers 1 and 2
        subplot(2,1,1)
    plot(time,zL1Tpos_ub,'b.',time,zL1Tpos_lb,'bx');
    grid on; title('Thawed Layer 1')
        subplot(2,1,2)
    plot(time,zL2Tpos_ub,'g.',time,zL2Tpos_lb,'gx');
    grid on; title('Thawed Layer 2')


    figure(2) % plot Thaw Layers 1,2,3
        subplot(3,1,1)
    plot(time,zL1Tpos_ub,'b.',time,zL1Tpos_lb,'bx');
    grid on; title('Thawed Layer 1')
        subplot(3,1,2)
    plot(time,zL2Tpos_ub,'g.',time,zL2Tpos_lb,'gx');
    grid on; title('Thawed Layer 2')
        subplot(3,1,3)
    plot(time,zL3Tpos_ub,'r.',time,zL3Tpos_lb,'rx');
    grid on; title('Thawed Layer 3')


    figure(3)
    hold on
    h1=plot(time,zL1Tpos_ub,'b.',time,zL1Tpos_lb,'bx');
    h2=plot(time,zL2Tpos_ub,'g.',time,zL2Tpos_lb,'gx');
    h3=plot(time,zL3Tpos_ub,'r.',time,zL3Tpos_lb,'rx');
    grid on; title('Thawed Layers'); legend([h1(1),h2(1),h3(1)],'Thaw Layer 1','Thaw Layer 2','Thaw Layer 3')
    hold off
end