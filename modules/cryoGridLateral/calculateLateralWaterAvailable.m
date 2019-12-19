function [ final_mat_worker, boundary_water_fluxes ] = calculateLateralWaterAvailable( PARA,GRID, wc, water_fluxes,boundary_water_fluxes )
% Function that calculates the real lateral fluxes based on the Darcy
% fluxes and the available water contents of the workers. The idea of the
% fonction is to use the "losing point of view" which considers that
% treating all the losing fluxes (water loss for the worker)  allows to
% treat all the fluxes between workers. Indeed,any gaining flux is a
% loosing flux from the reverse point of view. Therefore only dealing with
% the losing point of view should be enough.

final_mat_worker=zeros(numlabs,numlabs);

workerFluxes=water_fluxes(labindex,:);
lostWater=abs(sum(workerFluxes(workerFluxes<0)));
lostWater=lostWater + abs(min(boundary_water_fluxes,0));

if lostWater>0 % worker is loosing water
    
    % available water
    bottomBucketcTIndex = PARA.location.bottomBucketSoilcTIndex;
    soilType=GRID.soil.cT_soilType(1:bottomBucketcTIndex);
    fieldCapacity = zeros(size(soilType));
    for i=1:size(PARA.soil.soilTypes,1)
        fieldCapacity(soilType==i) = PARA.soil.soilTypes( i, 2 );
    end

    K_deltaSoil=GRID.general.K_delta(GRID.soil.cT_domain);   
    availableWaterVect=( wc(1:bottomBucketcTIndex)- fieldCapacity) .* K_deltaSoil(1:bottomBucketcTIndex);
    availableWaterVect( availableWaterVect < 0) = 0;
    availableWater=sum(availableWaterVect);
    fprintf('\t\t\tAvailable water (w%1.0f) : %3.2e m\n', labindex, availableWater)
    loss=find(workerFluxes<0);

    % Use the losing point of view
    if availableWater >= lostWater % when there is enough water, use the darcy fluxes
       scal_fact=1;
    else % when there is not enough water, scale lost fluxes so that it matches the available water
       scal_fact=availableWater/lostWater;
    end

    % Fill the worker matrix from the losing point of view
       for i=1:length(loss)
           final_mat_worker(labindex,loss(i))=scal_fact*workerFluxes(loss(i));
           final_mat_worker(loss(i),labindex)=scal_fact*(-1)*workerFluxes(loss(i))*PARA.ensemble.area(labindex)/PARA.ensemble.area(loss(i)) ...
               .* PARA.ensemble.hydraulic_contact_length(loss(i),labindex) ./ PARA.ensemble.hydraulic_contact_length(labindex, loss(i));
       end

    % Apply also to boundary fluxes
    boundary_water_fluxes=scal_fact*boundary_water_fluxes;
       
end

end



