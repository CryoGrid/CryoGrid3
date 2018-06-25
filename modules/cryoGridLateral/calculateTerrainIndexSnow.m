function terrain_index_snow_final = calculateTerrainIndexSnow(altitudes, weight)
%tsvd   
weight_fix = [1, 1]; % temporal fix  zzz ttt
weight=weight_fix;

snowCellSize = 0.025;
allAltitudes=[];
for i=1:size(weight,2)
   %allAltitudes=[allAltitudes repmat( altitudes(1,i) ,1, weight(1,i))];
   allAltitudes=[allAltitudes repmat( round(altitudes(1,i)./snowCellSize).*snowCellSize ,1, weight(1,i))];

end

terrain_index_snow=(allAltitudes-mean(allAltitudes))./std(allAltitudes);
terrain_index_snow=terrain_index_snow./sum(terrain_index_snow(terrain_index_snow>0)); %this must be normalized to conserve energy, check later with PARA.location.aerial_fraction
terrain_index_snow(isnan(terrain_index_snow))=0;

terrain_index_snow_final=[];
j=1;
for i=1:size(weight,2)
    terrain_index_snow_final=[terrain_index_snow_final terrain_index_snow(1,j)];
    j=j+weight(1,i);
end

%terrain_index_snow_final = terrain_index_snow_final./sum(terrain_index_snow_final(terrain_index_snow_final>0));