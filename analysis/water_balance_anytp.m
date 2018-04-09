clear
clc

% Script used to invistigate the goodness of the water conservation in mid
% March 2018

run='March14_5w5yr_snow_real1'; % Name of the global run (realisation number is the removed)
worker=4; % worker you want to work with <--------------------------------

% Load data
run=[run(1:end-1) num2str(worker)];
load([run '\' run '_output2016.mat'])

start='01-Aug-2015';
[~,start]=min(abs(OUT.timestamp-datenum(start)));
stop='01-Aug-2016';
[~,stop]=min(abs(OUT.timestamp-datenum(stop)));

% Initial and final states (m converted in mm with *1000), using
% parenthesis when you are not a cock sucker
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

lateral_snow    = 1000*nansum(OUT.WB.dr_lateralSnow(start:stop)); % WATCH OUT !!!!! if the x1000 is needed or not (run performed before the correction needs it, latter runs don't)
lateral_tot     = lateral_water+lateral_snow;

snow_excess     = nansum(OUT.WB.dr_excessSnow(start:stop)); % stored negtively
snow_melted     = nansum(OUT.WB.dr_snowmelt(start:stop)); % stored negtively

rain_surfRunoff  = nansum(OUT.WB.dr_surface(start:stop)); % stored negtively
external_flux   = nansum(OUT.WB.dr_external(start:stop)); % stored negtively

subli_condens   = nansum(OUT.WB.ds(start:stop)); % And condensation
evapotr         = nansum(OUT.WB.de(start:stop));

rain_lost = nansum(OUT.WB.dr_rain(start:stop));
lacking_water = nansum(OUT.WB.dm_lacking(start:stop));

input_output_sum         = inputs + lateral_tot + snow_excess + snow_melted + rain_surfRunoff + external_flux +  evapotr + rain_lost + subli_condens;

% Comparing with final-initial difference
delta_final_initial = final_water - initial_water;
balance_anomaly = input_output_sum - delta_final_initial;

% Calculating the brut worker fluxes to track back lateral exces or
% boundary water
lateral_waterFluxes_mat=OUT.lateral.water_fluxes(worker,:,(start:stop));
lateral_waterFluxes_mat=sum(lateral_waterFluxes_mat,2);
lateral_waterFluxes_mat=sum(lateral_waterFluxes_mat,3);

fprintf('Done\n')

%----------- Plotting ----------------------------------------------------

% Figure 1 - cumulative sum of the variables  
figure(1)
time=OUT.timestamp(start:stop);
plot(time,cumsum(OUT.WB.dp_rain(start:stop)+OUT.WB.dp_snow(start:stop)))

hold on
legend_txt{1,1}='Input';

plot(time,cumsum(OUT.WB.dr_lateral(start:stop)))
legend_txt{2,1}='Lateral_water';

plot(time,cumsum(OUT.WB.dr_snowmelt(start:stop)))
legend_txt{3,1}='snowmelt';

plot(time,cumsum(OUT.WB.ds(start:stop)))
legend_txt{4,1}='SubliCondens';

plot(time,cumsum(OUT.WB.de(start:stop)))
legend_txt{5,1}='ET';

plot(time,cumsum(OUT.WB.dr_rain(start:stop)))
legend_txt{6,1}='Lost rain';

plot(time,cumsum(OUT.WB.dr_DarcyReservoir(start:stop)))
legend_txt{7,1}='DarcyReservoir';

% plot(time,cumsum(OUT.WB.dr_surface))
% legend_txt{8,1}='surfaceRunoff';

% plot(time,cumsum(OUT.WB.dr_lateralSnow))
% legend_txt{9,1}='LateralSnow';

legend(legend_txt)
datetick
hold off

% Figure 2 - water table 
figure(2)
plot(time, OUT.location.water_table_altitude(start:stop))
hold on
datetick
hold off


% Figure 3 - brut left and right fluxes
figure(3)
if (worker>1 && worker<5)
    lateral_water_FluxesR=OUT.lateral.water_fluxes(worker,worker+1,(start:stop));
    lateral_water_FluxesR=permute(lateral_water_FluxesR, [3 2 1]);
    lateral_water_FluxesL=OUT.lateral.water_fluxes(worker,worker-1,(start:stop));
    lateral_water_FluxesL=permute(lateral_water_FluxesL, [3 2 1]);
    plot(time,cumsum(lateral_water_FluxesR))
    hold on
    plot(time,cumsum(lateral_water_FluxesL))
    legend_txt2{1,1}='Flux from uphill';
    legend_txt2{2,1}='Flux towards downhill';
    legend(legend_txt2)
    datetick
    hold off
end

if worker==1;
    lateral_water_FluxesR=OUT.lateral.water_fluxes(worker,worker+1,(start:stop));
    lateral_water_FluxesR=permute(lateral_water_FluxesR, [3 2 1]);
    plot(time,cumsum(lateral_water_FluxesR))
    hold on
    legend_txt2{1,1}='Flux from up hill';
    legend(legend_txt2)
    datetick
    hold off
end

if worker==5;
    lateral_water_FluxesL=OUT.lateral.water_fluxes(worker,worker-1,(start:stop));
    lateral_water_FluxesL=permute(lateral_water_FluxesL, [3 2 1]);
    plot(time,cumsum(lateral_water_FluxesL))
    hold on
    legend_txt2{1,1}='Flux towards down hill';
    legend(legend_txt2)
    datetick
    hold off
end

figure(4)
Val=[inputs lateral_tot snow_excess snow_melted rain_surfRunoff external_flux evapotr rain_lost subli_condens delta_final_initial];
Names={'inputs', 'lat tot', 'snow Xs', 'snow melt', 'rain surfRoff', 'ext flux', 'evapotr', 'rain lost', 'sub cond', 'd final init'};

[Val, I_sort]=sort(Val);
Names=Names(I_sort);

bar(1:length(Val), Val')
hold on
set(gca,'xticklabel',Names)
ylabel('balance contribution (mm)')
hold off
