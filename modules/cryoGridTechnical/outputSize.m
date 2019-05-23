function [ OUT ] = outputSize( savingVar, OUT  )
% Function that lower the amount of output of the code based on the
% variable PARA.technical.saving, inputed here as savingVar.
% savingVar adjusts the amount of files saved.
%    -1: Normal outputs
%     1: light outputs ONLY
%     2: light outputs, monthly mean values
%     3: light outputs, yearøy mean values
%     if higher than 1,2 or 3, it is to deal with the saving of the FINAL
%     structure and of the printing of the picture. So :
%     10, 20, 30 ~ 1, 2, 3 but with the picture printed 
%     100, 200, 300 ~ 1, 2, 3 but with the picture and the FINAL saved.

% Collect the good value
if savingVar~=-1;
    savingVar=sprintf('%1.2e',savingVar);
    savingVar=str2double(savingVar(1));
    assert(ismember(savingVar,[1 2 3]),'outputSize : Wrong savingVar')
    % Streamline
    OUT=rmfield(OUT,'debugging');
    OUT=rmfield(OUT,'EB');
    OUT=rmfield(OUT,'SEB');
    OUT=rmfield(OUT,'WB');
    OUT=rmfield(OUT,'lateral');
    OUT.snow=rmfield(OUT.snow,'outSnow_i');
    OUT.snow=rmfield(OUT.snow,'outSnow_a');
    OUT.snow=rmfield(OUT.snow,'outSnow_w');
    OUT.soil=rmfield(OUT.soil,'lakeFloor');
    OUT.soil=rmfield(OUT.soil,'soil');

    if savingVar==2; % Monthly means
        OUT=rmfield(OUT,'TIMESTEP');
        OUT.cryoGrid3=monthlymeanf_mat(OUT.timestamp, OUT.cryoGrid3);
        OUT.water=monthlymeanf_mat(OUT.timestamp, OUT.water);
        OUT.liquidWater=monthlymeanf_mat(OUT.timestamp, OUT.liquidWater);
        OUT.soil.topPosition=monthlymeanf_mat(OUT.timestamp, OUT.soil.topPosition);
        OUT.snow.topPosition=monthlymeanf_mat(OUT.timestamp, OUT.snow.topPosition);
        OUT.snow.botPosition=monthlymeanf_mat(OUT.timestamp, OUT.snow.botPosition);
        OUT.location.area=monthlymeanf_mat(OUT.timestamp, OUT.location.area);
        OUT.location.altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.altitude);
        OUT.location.soil_altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.soil_altitude);
        OUT.location.surface_altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.surface_altitude);
        OUT.location.infiltration_altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.infiltration_altitude);
        OUT.location.water_table_altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.water_table_altitude);
        OUT.location.infiltration_altitude_mean=monthlymeanf_mat(OUT.timestamp, OUT.location.infiltration_altitude_mean);
        OUT.location.water_table_altitude_mean=monthlymeanf_mat(OUT.timestamp, OUT.location.water_table_altitude_mean);
        OUT.location.pfTable_altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.pfTable_altitude);
        OUT.timestamp=monthlymeanf_mat(OUT.timestamp, OUT.timestamp);

    elseif savingVar==3 % Yearly means
        OUT=rmfield(OUT,'TIMESTEP');
        OUT.cryoGrid3=nanmean(OUT.cryoGrid3,2);
        OUT.water=nanmean(OUT.water,2);
        OUT.liquidWater=nanmean(OUT.liquidWater,2);
        OUT.soil.topPosition=nanmean(OUT.soil.topPosition);
        OUT.snow.topPosition=nanmean(OUT.snow.topPosition);
        OUT.snow.botPosition=nanmean(OUT.snow.botPosition);
        OUT.location.area=nanmean(OUT.location.area);
        OUT.location.altitude=nanmean(OUT.location.altitude);
        OUT.location.soil_altitude=nanmean(OUT.location.soil_altitude);
        OUT.location.surface_altitude=nanmean(OUT.location.surface_altitude);
        OUT.location.infiltration_altitude=nanmean(OUT.location.infiltration_altitude);
        OUT.location.water_table_altitude=nanmean(OUT.location.water_table_altitude);
        OUT.location.infiltration_altitude_mean=nanmean(OUT.location.infiltration_altitude_mean);
        OUT.location.water_table_altitude_mean=nanmean(OUT.location.water_table_altitude_mean);
        OUT.location.pfTable_altitude=nanmean(OUT.location.pfTable_altitude);
        OUT.timestamp=nanmean(OUT.timestamp);
    end
    
end