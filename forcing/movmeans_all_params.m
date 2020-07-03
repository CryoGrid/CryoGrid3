%clear all;

%% Load Forcing Data.
% NCAR = load('N:\permarisk\data\SNAP\MatFiles\NCAR_Forcing_1970-2100.mat');
% GFDL = load('N:\permarisk\data\SNAP\MatFiles\GFDL_Forcing_1970-2100.mat');
GFDL = load('GFDL_PrudhoeBay_RCP85_1970-2100');
NCAR = load('NCAR_PrudhoeBay_RCP85_1970-2100');
 
%% Calculate movmeans, masking the framing periods with incomplete input; saving data.
% % % timeFrame.name = {'decade', 'year', 'month'};
% % % timeFrame.duration = [87600, 8760, 730];
% % % parameters = {'Tair', 'rain', 'wind', 'snow', 'Sin', 'Lin', 'q', 'p'};
% % % models = {'NCAR', 'GFDL'};
% % % 
% % % for i = 1:numel(timeFrame.name)
% % %     for j = 1:numel(parameters)
% % %         for k = 1:numel(models)
% % %             pass
% % %         end
% % %     end
% % % end

% NCAR_decade.Tair = movmean(NCAR.FORCING.data.Tair, 87600);
% NCAR_decade.Tair(1:43800) = NaN;
% NCAR_decade.Tair(end-43800:end) = NaN;
% NCAR_decade.rain = movsum(NCAR.FORCING.data.rainfall, 87600)/240.0;
% NCAR_decade.rain(1:43800) = NaN;
% NCAR_decade.rain(end-43800:end) = NaN;
% NCAR_decade.wind = movmean(NCAR.FORCING.data.wind, 87600);
% NCAR_decade.wind(1:43800) = NaN;
% NCAR_decade.wind(end-43800:end) = NaN;
% NCAR_decade.snow = movsum(NCAR.FORCING.data.snowfall, 87600)/240.0;
% NCAR_decade.snow(1:43800) = NaN;
% NCAR_decade.snow(end-43800:end) = NaN;
% NCAR_decade.Sin = movmean(NCAR.FORCING.data.Sin, 87600);
% NCAR_decade.Sin(1:43800) = NaN;
% NCAR_decade.Sin(end-43800:end) = NaN;
% NCAR_decade.Lin = movmean(NCAR.FORCING.data.Lin, 87600);
% NCAR_decade.Lin(1:43800) = NaN;
% NCAR_decade.Lin(end-43800:end) = NaN;
% NCAR_decade.q = movmean(NCAR.FORCING.data.q, 87600);
% NCAR_decade.q(1:43800) = NaN;
% NCAR_decade.q(end-43800:end) = NaN;
% NCAR_decade.p = movmean(NCAR.FORCING.data.p, 87600);
% NCAR_decade.p(1:43800) = NaN;
% NCAR_decade.p(end-43800:end) = NaN;
% 
% GFDL_decade.Tair = movmean(GFDL.FORCING.data.Tair, 87600);
% GFDL_decade.Tair(1:43800) = NaN;
% GFDL_decade.Tair(end-43800:end) = NaN;
% GFDL_decade.rain = movsum(GFDL.FORCING.data.rainfall, 87600)/240.0;
% GFDL_decade.rain(1:43800) = NaN;
% GFDL_decade.rain(end-43800:end) = NaN;
% GFDL_decade.wind = movmean(GFDL.FORCING.data.wind, 87600);
% GFDL_decade.wind(1:43800) = NaN;
% GFDL_decade.wind(end-43800:end) = NaN;
% GFDL_decade.snow = movsum(GFDL.FORCING.data.snowfall, 87600)/240.0;
% GFDL_decade.snow(1:43800) = NaN;
% GFDL_decade.snow(end-43800:end) = NaN;
% GFDL_decade.Sin = movmean(GFDL.FORCING.data.Sin, 87600);
% GFDL_decade.Sin(1:43800) = NaN;
% GFDL_decade.Sin(end-43800:end) = NaN;
% GFDL_decade.Lin = movmean(GFDL.FORCING.data.Lin, 87600);
% GFDL_decade.Lin(1:43800) = NaN;
% GFDL_decade.Lin(end-43800:end) = NaN;
% GFDL_decade.q = movmean(GFDL.FORCING.data.q, 87600);
% GFDL_decade.q(1:43800) = NaN;
% GFDL_decade.q(end-43800:end) = NaN;
% GFDL_decade.p = movmean(GFDL.FORCING.data.p, 87600);
% GFDL_decade.p(1:43800) = NaN;
% GFDL_decade.p(end-43800:end) = NaN;
% 
% NCAR_year.Tair = movmean(NCAR.FORCING.data.Tair, 8760);
% NCAR_year.Tair(1:4380) = NaN;
% NCAR_year.Tair(end-4380:end) = NaN;
% NCAR_year.rain = movsum(NCAR.FORCING.data.rainfall, 8760)/24.0;
% NCAR_year.rain(1:4380) = NaN;
% NCAR_year.rain(end-4380:end) = NaN;
% NCAR_year.wind = movmean(NCAR.FORCING.data.wind, 8760);
% NCAR_year.wind(1:4380) = NaN;
% NCAR_year.wind(end-4380:end) = NaN;
% NCAR_year.snow = movsum(NCAR.FORCING.data.snowfall, 8760)/24.0;
% NCAR_year.snow(1:4380) = NaN;
% NCAR_year.snow(end-4380:end) = NaN;
% NCAR_year.Sin = movmean(NCAR.FORCING.data.Sin, 8760);
% NCAR_year.Sin(1:4380) = NaN;
% NCAR_year.Sin(end-4380:end) = NaN;
% NCAR_year.Lin = movmean(NCAR.FORCING.data.Lin, 8760);
% NCAR_year.Lin(1:4380) = NaN;
% NCAR_year.Lin(end-4380:end) = NaN;
% NCAR_year.q = movmean(NCAR.FORCING.data.q, 8760);
% NCAR_year.q(1:4380) = NaN;
% NCAR_year.q(end-4380:end) = NaN;
% NCAR_year.p = movmean(NCAR.FORCING.data.p, 8760);
% NCAR_year.p(1:4380) = NaN;
% NCAR_year.p(end-4380:end) = NaN;
% 
% GFDL_year.Tair = movmean(GFDL.FORCING.data.Tair, 8760);
% GFDL_year.Tair(1:4380) = NaN;
% GFDL_year.Tair(end-4380:end) = NaN;
% GFDL_year.rain = movsum(GFDL.FORCING.data.rainfall, 8760)/24.0;
% GFDL_year.rain(1:4380) = NaN;
% GFDL_year.rain(end-4380:end) = NaN;
% GFDL_year.wind = movmean(GFDL.FORCING.data.wind, 8760);
% GFDL_year.wind(1:4380) = NaN;
% GFDL_year.wind(end-4380:end) = NaN;
% GFDL_year.snow = movsum(GFDL.FORCING.data.snowfall, 8760)/24.0;
% GFDL_year.snow(1:4380) = NaN;
% GFDL_year.snow(end-4380:end) = NaN;
% GFDL_year.Sin = movmean(GFDL.FORCING.data.Sin, 8760);
% GFDL_year.Sin(1:4380) = NaN;
% GFDL_year.Sin(end-4380:end) = NaN;
% GFDL_year.Lin = movmean(GFDL.FORCING.data.Lin, 8760);
% GFDL_year.Lin(1:4380) = NaN;
% GFDL_year.Lin(end-4380:end) = NaN;
% GFDL_year.q = movmean(GFDL.FORCING.data.q, 8760);
% GFDL_year.q(1:4380) = NaN;
% GFDL_year.q(end-43800:end) = NaN;
% GFDL_year.p = movmean(GFDL.FORCING.data.p, 8760);
% GFDL_year.p(1:4380) = NaN;
% GFDL_year.p(end-4380:end) = NaN;
% 
% 
% NCAR_month.Tair = movmean(NCAR.FORCING.data.Tair, 730);
% NCAR_month.Tair(1:365) = NaN;
% NCAR_month.Tair(end-365:end) = NaN;
% NCAR_month.rain = movsum(NCAR.FORCING.data.rainfall, 730)/24.0; % mm/h
% NCAR_month.rain(1:365) = NaN;
% NCAR_month.rain(end-365:end) = NaN;
% NCAR_month.wind = movmean(NCAR.FORCING.data.wind, 730);
% NCAR_month.wind(1:365) = NaN;
% NCAR_month.wind(end-365:end) = NaN;
% NCAR_month.snow = movsum(NCAR.FORCING.data.snowfall, 730)/24.0;
% NCAR_month.snow(1:365) = NaN;
% NCAR_month.snow(end-365:end) = NaN;
% NCAR_month.Sin = movmean(NCAR.FORCING.data.Sin, 730);
% NCAR_month.Sin(1:365) = NaN;
% NCAR_month.Sin(end-365:end) = NaN;
% NCAR_month.Lin = movmean(NCAR.FORCING.data.Lin, 730);
% NCAR_month.Lin(1:365) = NaN;
% NCAR_month.Lin(end-365:end) = NaN;
% NCAR_month.q = movmean(NCAR.FORCING.data.q, 730);
% NCAR_month.q(1:365) = NaN;
% NCAR_month.q(end-365:end) = NaN;
% NCAR_month.p = movmean(NCAR.FORCING.data.p, 730);
% NCAR_month.p(1:365) = NaN;
% NCAR_month.p(end-365:end) = NaN;
% 
% GFDL_month.Tair = movmean(GFDL.FORCING.data.Tair, 730);
% GFDL_month.Tair(1:365) = NaN;
% GFDL_month.Tair(end-365:end) = NaN;
% GFDL_month.rain = movsum(GFDL.FORCING.data.rainfall, 730)/24.0;
% GFDL_month.rain(1:365) = NaN;
% GFDL_month.rain(end-365:end) = NaN;
% GFDL_month.wind = movmean(GFDL.FORCING.data.wind, 730);
% GFDL_month.wind(1:365) = NaN;
% GFDL_month.wind(end-365:end) = NaN;
% GFDL_month.snow = movsum(GFDL.FORCING.data.snowfall, 730)/24.0;
% GFDL_month.snow(1:365) = NaN;
% GFDL_month.snow(end-365:end) = NaN;
% GFDL_month.Sin = movmean(GFDL.FORCING.data.Sin, 730);
% GFDL_month.Sin(1:365) = NaN;
% GFDL_month.Sin(end-365:end) = NaN;
% GFDL_month.Lin = movmean(GFDL.FORCING.data.Lin, 730);
% GFDL_month.Lin(1:365) = NaN;
% GFDL_month.Lin(end-365:end) = NaN;
% GFDL_month.q = movmean(GFDL.FORCING.data.q, 730);
% GFDL_month.q(1:365) = NaN;
% GFDL_month.q(end-365:end) = NaN;
% GFDL_month.p = movmean(GFDL.FORCING.data.p, 730);
% GFDL_month.p(1:365) = NaN;
% GFDL_month.p(end-365:end) = NaN;
%
% output_path = 'N:\permarisk\data\SNAP\MatFiles\NCAR_decade.mat';
% save(output_path,'NCAR_decade');
% output_path = 'N:\permarisk\data\SNAP\MatFiles\GFDL_decade.mat';
% save(output_path,'GFDL_decade');
% output_path = 'N:\permarisk\data\SNAP\MatFiles\NCAR_year.mat';
% save(output_path,'NCAR_year');
% output_path = 'N:\permarisk\data\SNAP\MatFiles\GFDL_year.mat';
% save(output_path,'GFDL_year');
% output_path = 'N:\permarisk\data\SNAP\MatFiles\NCAR_month.mat';
% save(output_path,'NCAR_month');
% output_path = 'N:\permarisk\data\SNAP\MatFiles\GFDL_month.mat';
% save(output_path,'GFDL_month');

%%Load (averaged) Data.
% NCAR = load('N:\permarisk\data\SNAP\MatFiles\NCAR_Forcing_Fairbanks_1970-2100.mat');
% GFDL = load('N:\permarisk\data\SNAP\MatFiles\GFDL_Forcing_Fairbanks_1970-2100.mat');
NCAR = load('N:\permarisk\data\SNAP\MatFiles\NCAR_Forcing_Prudhoe_1970-2100.mat');
GFDL = load('N:\permarisk\data\SNAP\MatFiles\GFDL_Forcing_Prudhoe_1970-2100.mat');
load('N:\permarisk\data\SNAP\MatFiles\NCAR_decade.mat');
load('N:\permarisk\data\SNAP\MatFiles\GFDL_decade.mat');
load('N:\permarisk\data\SNAP\MatFiles\NCAR_year.mat');
load('N:\permarisk\data\SNAP\MatFiles\GFDL_year.mat');
load('N:\permarisk\data\SNAP\MatFiles\NCAR_month.mat');
load('N:\permarisk\data\SNAP\MatFiles\GFDL_month.mat');

%%Plot decade averaged data.
subplot(4,2,1);
plot(NCAR.FORCING.data.time, NCAR_decade.Tair, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.Tair, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('T [?C]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of temperatures');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,2);
plot(NCAR.FORCING.data.time, NCAR_decade.rain, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.rain, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('Rain [mm/h]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of rainfall');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,3);
plot(NCAR.FORCING.data.time, NCAR_decade.wind, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.wind, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('v [m/s]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of windspeeds');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,4);
plot(NCAR.FORCING.data.time, NCAR_decade.snow, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.snow, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('SWE [mm/h]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of snowfall');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,5);
plot(NCAR.FORCING.data.time, NCAR_decade.Sin, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.Sin, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('Incoming shortwave radiation [W/m2]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of downward shortwave flux');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,6);
plot(NCAR.FORCING.data.time, NCAR_decade.p, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.p, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('p [Pa]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of surface pressure');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,7);
plot(NCAR.FORCING.data.time, NCAR_decade.Lin, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.Lin, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('Incoming longwave radiation [W/m?]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of downward longwave flux');
grid(gca,'minor');
grid on;
hold off;

subplot(4,2,8);
plot(NCAR.FORCING.data.time, NCAR_decade.q, 'LineWidth', 1.5);
hold on;
plot(NCAR.FORCING.data.time, GFDL_decade.q, 'LineWidth', 1.5);
legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
ylabel('Specific humidity [%]');
xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
title('10-year-Moving mean of specific humidity');
grid(gca,'minor');
grid on;
hold off;




%%Plot yearly averaged data.
% subplot(4,2,1);
% plot(NCAR.FORCING.data.time, NCAR_year.Tair);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.Tair);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('T [?C]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of temperatures');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,2);
% plot(NCAR.FORCING.data.time, NCAR_year.rain);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.rain);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Rain [mm/h]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of rainfall');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,3);
% plot(NCAR.FORCING.data.time, NCAR_year.wind);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.wind);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('v [m/s]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of windspeeds');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,4);
% plot(NCAR.FORCING.data.time, NCAR_year.snow);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.snow);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('SWE [mm/h]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of snowfall');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,5);
% plot(NCAR.FORCING.data.time, NCAR_year.Sin);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.Sin);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Incoming shortwave radiation [W/m2]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of downward shortwave flux');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,6);
% plot(NCAR.FORCING.data.time, NCAR_year.p);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.p);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('p [Pa]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of surface pressure');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,7);
% plot(NCAR.FORCING.data.time, NCAR_year.Lin);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.Lin);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Incoming longwave radiation [W/m?]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of downward longwave flux');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,8);
% plot(NCAR.FORCING.data.time, NCAR_year.q);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_year.q);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Specific humidity [%]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-year-Moving mean of specific humidity');
% grid(gca,'minor');
% grid on;
% hold off;


%%Plot monthly averaged data.
% subplot(4,2,1);
% plot(NCAR.FORCING.data.time, NCAR_month.Tair);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.Tair);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('T [?C]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of temperatures');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,2);
% plot(NCAR.FORCING.data.time, NCAR_month.rain);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.rain);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Rain [mm/h]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of rainfall');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,3);
% plot(NCAR.FORCING.data.time, NCAR_month.wind);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.wind);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('v [m/s]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of windspeeds');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,4);
% plot(NCAR.FORCING.data.time, NCAR_month.snow);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.snow);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('SWE [mm/h]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of snowfall');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,5);
% plot(NCAR.FORCING.data.time, NCAR_month.Sin);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.Sin);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Incoming shortwave radiation [W/m2]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of downward shortwave flux');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,6);
% plot(NCAR.FORCING.data.time, NCAR_month.p);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.p);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('p [Pa]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of surface pressure');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,7);
% plot(NCAR.FORCING.data.time, NCAR_month.Lin);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.Lin);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Incoming longwave radiation [W/m?]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of downward longwave flux');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,8);
% plot(NCAR.FORCING.data.time, NCAR_month.q);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL_month.q);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Specific humidity [%]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('1-month-Moving mean of specific humidity');
% grid(gca,'minor');
% grid on;
% hold off;



%%Plot original data.
% subplot(4,2,1);
% plot(NCAR.FORCING.data.time, NCAR.FORCING.data.Tair);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL.FORCING.data.Tair);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('T [?C]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly temperatures');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,2);
% bar(NCAR.FORCING.data.time, [NCAR.FORCING.data.rainfall/24.0 GFDL.FORCING.data.rainfall/24.0]);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Rain [mm/h]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly rainfall');
% grid(gca,'minor');
% grid on;
% 
% subplot(4,2,3);
% plot(NCAR.FORCING.data.time, NCAR.FORCING.data.wind);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL.FORCING.data.wind);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('v [m/s]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly windspeeds');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,4);
% bar(NCAR.FORCING.data.time, [NCAR.FORCING.data.snowfall/24.0 GFDL.FORCING.data.snowfall/24.0]);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('SWE [mm/h]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly snowfall');
% grid(gca,'minor');
% grid on;
% 
% subplot(4,2,5);
% plot(NCAR.FORCING.data.time, NCAR.FORCING.data.Sin);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL.FORCING.data.Sin);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Incoming shortwave radiation [W/m2]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly downward shortwave flux');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,6);
% plot(NCAR.FORCING.data.time, NCAR.FORCING.data.p);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL.FORCING.data.p);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('p [Pa]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly surface pressure');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,7);
% plot(NCAR.FORCING.data.time, NCAR.FORCING.data.Lin);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL.FORCING.data.Lin);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Incoming longwave radiation [W/m?]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly downward longwave flux');
% grid(gca,'minor');
% grid on;
% hold off;
% 
% subplot(4,2,8);
% plot(NCAR.FORCING.data.time, NCAR.FORCING.data.q);
% hold on;
% plot(NCAR.FORCING.data.time, GFDL.FORCING.data.q);
% legend('CCSM 4.0 (\it{NCAR})', 'CM 3.0 (\it{GFDL})');
% ylabel('Specific humidity [%]');
% xlim([NCAR.FORCING.data.time(1) NCAR.FORCING.data.time(end)]);
% title('Hourly specific humidity');
% grid(gca,'minor');
% grid on;
% hold off;