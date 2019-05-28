function [ index, OUT ] = pfTable( OUT, mode )
% Find the permafrost table from a cryoGrid3 matrix

cryoGrid3=OUT.cryoGrid3;
[~,nbt]=size(cryoGrid3);

frozen = cryoGrid3<=0;
frozen2 = sum(frozen,2);
frozen3 = frozen2==nbt;
frozen4=find(frozen3==1);

if isempty(frozen4)
    index=NaN;
else
    index=frozen4(1);
end

if nargin==1;
    mode=0;
end

if mode==1 && isnan(index)
    nval=size(OUT.cryoGrid3,2);
    OUT.location.pfTable_altitude=nan(nval,1);
end

end