clear
close all

% Script used to invistigate the goodness of the water conservation in mid
% March 2018
% Amended in April because now the lateral excess water is driven toward
% the forming of a small pond, so now enter into the budget calculation,
% and excess appears in runoff.

run='April30_5yr_fig'; % Name of the global run (realisation number is the removed)
addpath \\lagringshotell\geofag\projects\coup\Leo\toolbox
name_str = runFiles( run );

worker=1; % worker you want to work with <--------------------------------

% Load data
load(fullfile(run, name_str(worker).output{end}));

start_str='01-Aug-2015';
[~,start]=min(abs(OUT.timestamp-datenum(start_str)));
stop_str='01-Aug-2016';
[~,stop]=min(abs(OUT.timestamp-datenum(stop_str)));

StoreAnomaly=zeros(stop-start,1);
StoreDeltaFI=zeros(stop-start,1);
StoreInOutSum=zeros(stop-start,1);

for i_stoptoplus=1:(stop-start+1)-1;
    i_stop=start+i_stoptoplus;
    
    % Initial and final states (m converted in mm with *1000), using
    % parenthesis when you are not a cock sucker
    K_delta_i=(-OUT.debugging.K_grid(1:end-1,start) + OUT.debugging.K_grid(2:end,start));
    initial_water=(nansum(K_delta_i.*OUT.water(:,start)) + nansum(OUT.snow.outSnow_i(:,start)) + nansum(OUT.snow.outSnow_w(:,start)))*1000 + OUT.WB.water2pool(start); %*1000; % *1000 to have mm check if needed !!
    K_delta_f=(-OUT.debugging.K_grid(1:end-1,i_stop) + OUT.debugging.K_grid(2:end,i_stop));
    final_water=(nansum(K_delta_f.*OUT.water(:,i_stop)) + nansum(OUT.snow.outSnow_i(:,i_stop)) + nansum(OUT.snow.outSnow_w(:,i_stop)))*1000 + OUT.WB.water2pool(i_stop); % *1000; % *1000 to have mm check if needed !!
    clear K_delta_i K_delta_f
    
    % Inputs : snow and rain (mm)
    inputs= nansum(OUT.WB.dp_rain(start:i_stop))+nansum(OUT.WB.dp_snow(start:i_stop));
    
    % Variations (all stored in mm)
    lateral_water = nansum(OUT.WB.dr_lateralWater(start:i_stop));
    lateral_darcy   = nansum(OUT.WB.dr_DarcyReservoir(start:i_stop));
    lateral_excessWater  = nansum(OUT.WB.dr_lateralExcess(start:i_stop));
    
    lateral_snow    = nansum(OUT.WB.dr_lateralSnow(start:i_stop)); % *1000; % WATCH OUT !!!!! if the x1000 is needed or not (run performed before the correction needs it, latter runs don't)
    lateral_tot     = lateral_excessWater + lateral_water+lateral_snow;
    
    snow_excess     = nansum(OUT.WB.dr_excessSnow(start:i_stop)); % stored negtively
    snow_melted     = 0; % nansum(OUT.WB.dr_snowmelt(start:i_stop)); % stored negtively, reinvested in GRID.residualWater, and finally in surface runoff, not to be counted twice
    
    rain_surfRunoff  = nansum(OUT.WB.dr_surface(start:i_stop)); % stored negtively
    external_flux   = nansum(OUT.WB.dr_external(start:i_stop)); % stored negtively
    
    subli_condens   = nansum(OUT.WB.ds(start:i_stop)); % And condensation
    evapotr         = nansum(OUT.WB.de(start:i_stop));
    
    rain_lost = 0; % nansum(OUT.WB.dr_rain(start:i_stop)); % reinvested in GRID.residualWater, and finally in surface runoff, not to be counted twice
    lacking_water = nansum(OUT.WB.dm_lacking(start:i_stop));
    
    input_output_sum         =  inputs + lateral_tot + snow_excess + snow_melted + rain_surfRunoff + external_flux +  evapotr + rain_lost + subli_condens;
    
    % Comparing with final-initial difference
    delta_final_initial = final_water - initial_water;
    balance_anomaly = input_output_sum - delta_final_initial;
    
    StoreAnomaly(i_stoptoplus,1)=balance_anomaly;
    StoreDeltaFI(i_stoptoplus,1)=delta_final_initial;
    StoreInOutSum(i_stoptoplus,1)=input_output_sum;
    
end
fprintf('run : %s\tStart : %s\t Stop : %s\tworker : %1.0f\n',run, start_str, stop_str,worker)
fprintf('Balance anomaly     : %2.1f mm\n', StoreAnomaly(end))
fprintf('Delta final-initial : %2.1f mm\n', delta_final_initial)
fprintf('Sum in and outputs  : %2.1f mm\n', input_output_sum)

figure(1)
figure('Position', [40, 630, 800,400]);
plot(OUT.timestamp(start+1:stop),StoreAnomaly)
hold on
plot(OUT.timestamp(start+1:stop),StoreDeltaFI) %-124.7)
plot(OUT.timestamp(start+1:stop),StoreInOutSum) %-124.7)
% plot(OUT.timestamp(start+1:stop),OUT.snow.topPosition(start+1:stop).*1000)
% plot(OUT.timestamp(start+1:stop),(OUT.cryoGrid3(83,start+1:stop)>0).*10,'.')
% plot(OUT.timestamp(start+1:stop),OUT.WB.residualWater(start+1:stop).*1000)
xlabel('time')
ylabel('water budget (mm)')
legend('Anomaly','Delta FI','SumInOut') %,'snow', 'top soil unfrozen','residualWater')
datetick
hold off

fprintf('Done\n')

%----------- Plotting ----------------------------------------------------

% Figure 1 - cumulative sum of the variables  
figure(2)
figure('Position', [990, 630, 800,400]);
time=OUT.timestamp(start:stop);
plot(time,cumsum(OUT.WB.dp_rain(start:stop)))

hold on
legend_txt{1,1}='rain input';

plot(time,cumsum(OUT.WB.dp_snow(start:stop)))
legend_txt{2,1}='snow_input';

plot(time,(cumsum(OUT.WB.dr_lateralSnow(start:i_stop))+cumsum(OUT.WB.dr_lateralExcess(start:stop))+cumsum(OUT.WB.dr_lateralWater(start:stop)))) % see if *1000 needed for snow
legend_txt{3,1}='lateral_tot';

plot(time,cumsum(OUT.WB.dr_excessSnow(start:stop)),'.')
legend_txt{4,1}='Xs snow';

plot(time,cumsum(OUT.WB.dr_surface(start:stop)),'.')
legend_txt{5,1}='surfaceRunoff';

plot(time,cumsum(OUT.WB.de(start:stop)))
legend_txt{6,1}='ET';

plot(time,cumsum(OUT.WB.dr_rain(start:stop)))
legend_txt{7,1}='lost rain';

plot(time,cumsum(OUT.WB.ds(start:stop)))
legend_txt{8,1}='SubliCondens';

% plot(time,cumsum(OUT.WB.dr_external(start:stop)))
% legend_txt{9,1}='external';
% plot(time,cumsum(OUT.WB.dr_snowmelt(start:stop)))
% legend_txt{10,1}='snowmelt';
% plot(time,cumsum(OUT.WB.dr_lateralSnow))
% legend_txt{11,1}='LateralSnow';
% plot(time,cumsum(OUT.WB.dr_DarcyReservoir(start:stop)))
% legend_txt{12,1}='DarcyReservoir';
% plot(time,cumsum(OUT.WB.dr_lateralExcess(start:stop)),'.')
% legend_txt{13,1}='Lateral xs water';

legend(legend_txt)
datetick
hold off

% Figure 2 - water table 
figure(3)
figure('Position', [40, 110, 800,400]);
plot(time, OUT.location.water_table_altitude(start:stop))
hold on
datetick
hold off


% Figure 3 - brut left and right fluxes
% figure(3)
% if (worker>1 && worker<5)
%     lateral_water_FluxesR=OUT.lateral.water_fluxes(worker,worker+1,(start:stop));
%     lateral_water_FluxesR=permute(lateral_water_FluxesR, [3 2 1]);
%     lateral_water_FluxesL=OUT.lateral.water_fluxes(worker,worker-1,(start:stop));
%     lateral_water_FluxesL=permute(lateral_water_FluxesL, [3 2 1]);
%     plot(time,cumsum(lateral_water_FluxesR))
%     hold on
%     plot(time,cumsum(lateral_water_FluxesL))
%     legend_txt2{1,1}='Flux from uphill';
%     legend_txt2{2,1}='Flux towards downhill';
%     legend(legend_txt2)
%     datetick
%     hold off
% end

% if worker==1;
%     lateral_water_FluxesR=OUT.lateral.water_fluxes(worker,worker+1,(start:stop));
%     lateral_water_FluxesR=permute(lateral_water_FluxesR, [3 2 1]);
%     plot(time,cumsum(lateral_water_FluxesR))
%     hold on
%     legend_txt2{1,1}='Flux from up hill';
%     legend(legend_txt2)
%     datetick
%     hold off
% end
% 
% if worker==2;
%     lateral_water_FluxesL=OUT.lateral.water_fluxes(worker,worker-1,(start:stop));
%     lateral_water_FluxesL=permute(lateral_water_FluxesL, [3 2 1]);
%     plot(time,cumsum(lateral_water_FluxesL))
%     hold on
%     legend_txt2{1,1}='Flux towards down hill';
%     legend(legend_txt2)
%     datetick
%     hold off
% end

figure(4)
figure('Position', [990, 110, 800,400]);
Val=[inputs lateral_tot snow_excess snow_melted rain_surfRunoff external_flux evapotr rain_lost subli_condens -delta_final_initial];
Names={'inputs', 'lat tot', 'snow Xs', 'snow melt', 'rain surfRoff', 'ext flux', 'evapotr', 'rain lost', 'sub cond', '-1* d final init'};

[Val, I_sort]=sort(Val);
Names=Names(I_sort);

bar(1:length(Val), Val')
hold on
set(gca,'xticklabel',Names)
ylabel('balance contribution (mm)')
hold off
