function [ OUT ] = outputSize( savingVar, OUT  )
% Function that lower the amount of output of the code based on the
% variable PARA.technical.saving, inputed here as savingVar.
% savingVar adjusts the amount of files saved.
%    -1: Normal outputs
%     1: light outputs ONLY
%     2: ultra light outputs ONLY
%    10: light outputs + Plot
%    20: ultra light outputs + Plot
%   100: light outputs + Plot + FINAL
%   200: ultra light outputs + Plot + FINAL

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
        OUT.location.pfTable_altitude=monthlymeanf_mat(OUT.timestamp, OUT.location.OUT.location.pfTable_altitude);

        OUT=rmfield(OUT,'timestamp');
        OUT=rmfield(OUT,'TIMESTEP');

    elseif savingVar==3 % Yearly means
        OUT=rmfield(OUT,'timestamp');
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
        OUT.location.pfTable_altitude=nanmean(OUT.location.OUT.location.pfTable_altitude);
    end
    
end