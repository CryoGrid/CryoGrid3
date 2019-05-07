function [ stratiOut ] = stratiXice( initialStratigraphy, palsaHeight, ActiveL  )
% Function that add excess ice in the stratigraphy given the height of the
% palsa. 

% Define normal mineral and organic proportions
propMin=initialStratigraphy(2,3) / ( initialStratigraphy(2,3) + initialStratigraphy(2,4) );
propOrg=1-propMin;

% New contents
soilDepth=initialStratigraphy(3,1);
solidTot = ( 1 - initialStratigraphy(2,6)) * ((soilDepth - ActiveL - palsaHeight) / (soilDepth - ActiveL));
newMin = solidTot * propMin;
newOrga = solidTot * propOrg;
newIcew=1-newMin-newOrga;

% Create output
stratiOut=[initialStratigraphy(1,:);...
           ActiveL newIcew newMin newOrga initialStratigraphy(2,5:6);...
           initialStratigraphy(3:end,:)];
end
