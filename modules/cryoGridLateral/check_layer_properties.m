function check_layer_properties(PARA)
% make consistency check for soil layer properties

num_layers = size(PARA.soil.layer_properties,1);               
%check_sum = sum( PARA.soil.layer_properties(:,3)+PARA.soil.layer_properties(:,4)+PARA.soil.layer_properties(:,6) );
%assert(check_sum == num_layers,[' min + org + porosity do not add up to 1 in all layers  for worker ',num2str(labindex)])
PARA = loadSoilTypes( PARA );     
wc=PARA.soil.layer_properties(:,2); ms=PARA.soil.layer_properties(:,3); org=PARA.soil.layer_properties(:,4); soilType = PARA.soil.layer_properties(:,5); por = PARA.soil.layer_properties(:,6);

for i=1:size(PARA.soil.soilTypes,1)
	FC(soilType==i) = PARA.soil.soilTypes( i, 2 ); % determine field capacity
end
% check if field Capacity smaller than porosity
for i=1:num_layers-1 % discard lowest bedrock layer for consistency check of field capacity
%     if( FC(i) >= PARA.soil.layer_properties(i,6) );  disp([' FC not smaller than porosity at layer ',num2str(i)]); end
    assert( FC(i) < por(i),[' FC not smaller than porosity at layer ',num2str(i)])
end

% check if layer constituents add up to one
for i=1:num_layers
    if( wc(i)<=por(i) )
        assert( ms(i)+org(i)+por(i) == 1,[' ms+org+porosity do not equal 1 for layer ',num2str(i)] )
    else % xice
        assert( wc(i)+ms(i)+org(i) == 1,[' XICE: wc+ms+org do not equal 1 for layer ',num2str(i)] )
    end
end

   
  
    
   
