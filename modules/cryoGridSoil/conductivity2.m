function conductivity2 = conductivity2(water, ice, mineral, organic, PARA)

ka = PARA.constants.k_a; %0.025;       %air [Hillel(1982)]
kw = PARA.constants.k_w; %0.57;        %water [Hillel(1982)]
ko = PARA.constants.k_o; %0.25;        %organic [Hillel(1982)]
km = PARA.constants.k_m; %soil.kh_bedrock;     %mineral 
ki = PARA.constants.k_i; %2.2;         %ice [Hillel(1982)]

air=1-water-ice-mineral-organic;

conductivity2= (water.* kw.^0.5 + ice.* ki.^0.5 + mineral.* km.^0.5 + organic.* ko.^0.5 + air.* ka.^0.5).^2;

