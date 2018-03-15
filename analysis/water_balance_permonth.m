clear
clc

% Script used to invistigate the goodness of the water conservation in mid
% March 2018

run='March14_5w5yr_full_real3';

% Prepare the calculation for all the workers and all the months
balance_tab=zeros(5,12);
months={'01-Aug-2015' '01-Sep-2015' '01-Oct-2015' '01-Nov-2015' '01-Dec-2015' '01-Jan-2016' '01-Feb-2016' '01-Mar-2016' '01-Apr-2016' '01-May-2016' '01-Jun-2016' '01-Jul-2016' '01-Aug-2016'};

for worker=1:5;
    
    run=[run(1:end-1) num2str(worker)];
    load([run '\' run '_output2016.mat'])
    
    for i_month=1:12;
        
        % Dates
        start=months{i_month};
        [~,start]=min(abs(OUT.timestamp-datenum(start)));
        stop=months{i_month+1};
        [~,stop]=min(abs(OUT.timestamp-datenum(stop)));
        
        % Initial and final states (m converted in mm with *1000)
        K_delta_i=(-OUT.debugging.K_grid(1:end-1,start) + OUT.debugging.K_grid(2:end,start));
        initial_water=(nansum(K_delta_i.*OUT.water(:,start)) + nansum(OUT.snow.outSnow_i(:,start)) + nansum(OUT.snow.outSnow_w(:,start)))*1000; % *1000 to have mm
        K_delta_f=(-OUT.debugging.K_grid(1:end-1,stop) + OUT.debugging.K_grid(2:end,stop));
        final_water=(nansum(K_delta_f.*OUT.water(:,stop)) + nansum(OUT.snow.outSnow_i(:,stop)) + nansum(OUT.snow.outSnow_w(:,stop)))*1000; % *1000 to have mm
        clear K_delta_i K_delta_f
        
        % Inputs : snow and rain (mm)
        inputs= nansum(OUT.WB.dp_rain(start:stop))+nansum(OUT.WB.dp_snow(start:stop));        
        
        % Variations (all stored in mm)
        lateral_water = nansum(OUT.WB.dr_lateral(start:stop));
        lateral_darcy   = nansum(OUT.WB.dr_DarcyReservoir(start:stop));
        lateral_excessWater  = nansum(OUT.WB.dr_lateralExcess(start:stop));
        
        lateral_snow    = 1000*nansum(OUT.WB.dr_lateralSnow(start:stop));
        lateral_tot     = lateral_water+lateral_snow;
        
        snow_excess     = nansum(OUT.WB.dr_excessSnow(start:stop)); % stored negtively
        snow_melted     = nansum(OUT.WB.dr_snowmelt(start:stop)); % stored negtively
        
        rain_surfRunoff = nansum(OUT.WB.dr_surface(start:stop)); % stored negtively
        external_flux   = nansum(OUT.WB.dr_external(start:stop)); % stored negtively
        
        subli_condens   = nansum(OUT.WB.ds(start:stop)); % And condensation
        evapotr         = nansum(OUT.WB.de(start:stop));
        
        rain_lost       = nansum(OUT.WB.dr_rain(start:stop));
        lacking_water   = nansum(OUT.WB.dm_lacking(start:stop));
        
        input_output_sum = inputs + lateral_tot + snow_excess + snow_melted + rain_surfRunoff + external_flux +  evapotr + rain_lost + subli_condens;
        delta_final_initial = final_water - initial_water;
        
        balance_anomaly = input_output_sum - delta_final_initial;
        
        balance_tab(worker,i_month)=balance_anomaly;
    end
    
end

fprintf('Worker %2g, annual balance anomaly : %3.2e mm\n',[[1:5] ; sum(balance_tab,2)'])

% clear worker start stop i i_month

%----- Plotting ----------------------------------------------------------

% Preparing the legend
for i=1:13;
    months{i}=months{i}(4:end);
    months{i}(end-4:end-2)=[];
end

% Plotting per month for 2 selected workers
bar(1:12, [balance_tab(5,:)' balance_tab(2,:)'])
hold on
set(gca,'xticklabel',months(1:end-1))
ylabel('balance anomaly (mm)')
legend({'Top of the slope' 'Down the slope'})
hold off