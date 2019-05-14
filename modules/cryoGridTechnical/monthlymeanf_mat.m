function [ monthlyval ] = monthlymeanf_mat( time,val,calcway)
% Function that uses the monthlymeanf on YEARLY outputs matrix to compute
% the monthy mean of each line. Input format is specific so check this out.
% (Example : time alone line and depth along column, ...)
% calcway==0 -> mean
% calcway==1 -> nanmean

if nargin<=2;
    calcway=1;
end

[nl,nc]=size(val);
month_index=1:12;
starting_month=month(time(1));
month_index=circshift(month_index,[0,-(starting_month-1)]);

if min(nl,nc)==1; % Vector
    monthlyval=nan(1,12);
    for i_m=1:12;
        time_li=month(time)==month_index(i_m);
        if calcway==0;
            monthlyval(i_m)=mean(val(time_li));
        else
            monthlyval(i_m)=nanmean(val(time_li));
        end
    end
else % Matrix
    monthlyval=nan(nl,12);
    % nl=size(val,1);
    for i_m=1:12;
        time_li=month(time)==month_index(i_m);
        time_li=logical(repmat(time_li',[nl,1]));
        nc=max(sum(time_li,2));
        extract=val(time_li);
        extract=reshape(extract,[nl nc]);
        if calcway==0;
            monthlyval(:,i_m)=mean(extract,2);
        else
            monthlyval(:,i_m)=nanmean(extract,2);
        end
    end
    
end

end