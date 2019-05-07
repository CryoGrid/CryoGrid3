function [ stratiOut ] = stratiXiceOld( initialStratigraphy, palsaHeight, ActiveL  )
% Function that add excess ice in the stratigraphy given the height of the
% palsa. 

% Consider how the varable are beeing taken
soilDepth=initialStratigraphy(3,1);
initialStratigraphy(2,1)=ActiveL;
expanded=soilDepth-ActiveL;

% New contents
newMin =initialStratigraphy(2,3)*expanded/(expanded+palsaHeight);
newOrga=initialStratigraphy(2,4)*expanded/(expanded+palsaHeight);
newIcew=1-newMin-newOrga;

% Create output
stratiOut=[initialStratigraphy(1,:);...
           ActiveL newIcew newMin newOrga initialStratigraphy(2,5:6);...
           initialStratigraphy(3:end,:)];
       
end

