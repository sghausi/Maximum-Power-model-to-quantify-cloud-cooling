%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%   To Implement the maximum power approach and quantify cloud cooling   %
%                         By Sarosh Alam Ghausi                          %
%             Max Planck Institute for Biogeochemistry, Jena - Germany   %     
%                   Email: sghausi@bgc-jena.mpg.de                       %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% **** Reading variables from NASA-CERES nc file  ****
lat_globe = ncread('CERES_globe.nc','lat');
lon_globe = ncread('CERES_globe.nc','lon');
lon_g2 = rem((lon_globe+180),360)-180; %converting long from (0:360) into (-180:180)

Net_sw_globe_all = ncread('CERES_globe.nc','sfc_net_sw_all_mon');
Net_sw_globe_clr = ncread('CERES_globe.nc','sfc_net_sw_clr_t_mon');
lwdown_globe_all = ncread('CERES_globe.nc','sfc_lw_down_all_mon');
lwdown_globe_clr = ncread('CERES_globe.nc','sfc_lw_down_clr_t_mon');
lwup_globe_all = ncread('CERES_globe.nc','sfc_lw_up_all_mon');
lwup_globe_clr = ncread('CERES_globe.nc','sfc_lw_up_clr_t_mon');
cld_area_globe = ncread('CERES_globe.nc','cldarea_total_daynight_mon');
solar_incom_globe = ncread('CERES_globe.nc','solar_mon');
sw_toa_globe_all = ncread('CERES_globe.nc','toa_sw_all_mon');
sw_toa_globe_clr = ncread('CERES_globe.nc','toa_sw_clr_t_mon');
lw_toa_globe_all = ncread('CERES_globe.nc','toa_lw_all_mon');
lw_toa_globe_clr = ncread('CERES_globe.nc','toa_lw_clr_t_mon');

sw_down_all = ncread('ceres_globe_sw.nc','sfc_sw_down_all_mon');
sw_up_all = ncread('ceres_globe_sw.nc','sfc_sw_up_all_mon');
sw_down_clr = ncread('ceres_globe_sw.nc','sfc_sw_down_clr_t_mon');
sw_up_clr = ncread('ceres_globe_sw.nc','sfc_sw_up_clr_t_mon');

% Download daily CERES daya from https://asdc.larc.nasa.gov/project/CERES/CER_SYN1deg-Day_Terra-Aqua-MODIS_Edition4A

%% ****** Processing data for each grid ******

load('globe_grid.mat') 
%load the global land mask generated using IGBP land-cover classification

for s = 1:length(grid_mask)

      globe_sw_all{s,1} = squeeze(Net_sw_globe_all(lon_index,lat_index,:));
      globe_sw_clr{s,1} = squeeze(Net_sw_globe_clr(lon_index,lat_index,:));
      globe_ldw_all{s,1} = squeeze(lwdown_globe_all(lon_index,lat_index,:));
      globe_ldw_clr{s,1} = squeeze(lwdown_globe_clr(lon_index,lat_index,:));
      globe_lwup_all{s,1} = squeeze(lwup_globe_all(lon_index,lat_index,:));
      globe_lwup_clr{s,1} = squeeze(lwup_globe_clr(lon_index,lat_index,:));
      globe_cld_area{s,1} = squeeze(cld_area_globe(lon_index,lat_index,:));
      globe_lwtoa_all{s,1} = squeeze(lw_toa_globe_all(lon_index,lat_index,:));
      globe_lwtoa_clr{s,1} = squeeze(lw_toa_globe_clr(lon_index,lat_index,:));
      globe_solarincom{s,1} = squeeze(solar_incom_globe(lon_index,lat_index,:));
      globe_sw_toa{s,1} = squeeze(sw_toa_globe_all(lon_index,lat_index,:));
      globe_sw_toa_clr{s,1} = squeeze(sw_toa_globe_clr(lon_index,lat_index,:));
      globe_sw_down_clr{s,1} = squeeze(sw_down_clr(lon_index,lat_index,:));
      globe_sw_down_all{s,1} = squeeze(sw_down_all(lon_index,lat_index,:));
      globe_sw_up_clr{s,1} = squeeze(sw_up_clr(lon_index,lat_index,:));
      globe_sw_up_all{s,1} = squeeze(sw_up_all(lon_index,lat_index,:));

      % This creates a structure array for each variable such that each
      % cell represent a grid and its time series
disp(s)
  end

 %% ***** Implementing the Maximum Power Model with heat storage for all-sky *****

for s = 1:length(grid_mask)
    globe_rnet{s,1} = globe_sw_all{s,1} + globe_ldw_all{s,1} - globe_lwup_all{s,1}; % defines net radiation
    globe_rlnet{s,1} = globe_lwup_all{s,1} - globe_ldw_all{s,1}; %defines net longwave radiation 
    du{s,1} = globe_sw_all{s,1} - globe_lwtoa_all{s,1}; % defines heat storage
    du_mean(s,1) = mean(du{s,1});
    disp(s)
end
 
 for s = 1:length(grid_mask)
Rin{s,1} = globe_sw_all{s,1} + globe_ldw_all{s,1}; % defines total radiative heating at the surface
globe_rlnet{s,1} = globe_lwup_all{s,1} - globe_ldw_all{s,1} ;
T_rad_all{s,1} = (globe_lwtoa_all{s,1}/(1* 5.67 * 10^-8)).^0.25; % defines the radiative temperature of atmosphere
T_rad_clr{s,1} = (globe_lwtoa_clr{s,1}/(1* 5.67 * 10^-8)).^0.25;
T_estimated{s,1} = (globe_lwup_all{s,1}/(5.67 * 10^-8)).^0.25; % Estimated skin temperature from CERES
end


for s = 1:length(grid_mask)
      for i = 1:length(days) % Running for each day 

jt = []; 
m = [];
idx = [];
G = []; 

if Rin{s,1}(i,1)  - globe_lwtoa_all{s,1}(i,1) > 0  % cases that does not undergo inversion 

jt = (1:(Rin{s,1}(i,1)  - globe_lwtoa_all{s,1}(i,1)))';

for w = 1:length(jt) % calculate the power and temperature for every possible value of J
ts(w,1) = ((Rin{s,1}(i,1)  - jt(w,1))/(5.67 * 10^-8)).^0.25; % From energy balance
G(w,1) = (jt(w,1)- du{s,1}(i,1)) * ((ts(w,1)/T_rad_all{s,1}(i,1)) - 1); % Carnot limit for dissipative heat engine
end

% Numerically solving for dG/dJ = 0

[m,idx] = max(G); % idx represent turbulent flux (J) at max power
power{s,1}(i,1) = m;
turb_all{s,1}(i,1) = idx;
%LE_mp{s,1}(i,1) = (slope{s,1}(i,1)/(slope{s,1}(i,1) + 65)) * (turb{s,1}(i,1));
%H_ceres{s,1}(i,1) = (65/(slope{s,1}(i,1) + 65)) * (turb{s,1}(i,1));
Tmp_all{s,1}(i,1) = ((Rin{s,1}(i,1) - turb_all{s,1}(i,1))/(5.67 * 10^-8)).^0.25;
jt = []; 
m = [];
idx = [];
G = [];
else % Removing the cases that undergoes inversion

idx = 0;
Tmp_all{s,1}(i,1) = ((Rin{s,1}(i,1) - globe_rnet{s,1}(i,1))/(5.67 * 10^-8)).^0.25;
%LE_mp{s,1}(i,1) = (slope{s,1}(i,1)/(slope{s,1}(i,1) + 65)) * (globe_rnet{s,1}(i,1));
turb_all{s,1}(i,1) = globe_rnet{s,1}(i,1);
jt = [];
m = [];
idx = [];
G = [];  
end
    end
    disp(s)
end

%% *** Implementing max power heat storage model for clear sky conditions ***

for s = 1:length(grid_mask)
Rin_clr{s,1} = globe_sw_clr{s,1} + globe_ldw_clr{s,1};
 du_clr{s,1} = globe_sw_clr{s,1} - globe_lwtoa_clr{s,1};
 globe_rnet_clr{s,1} = globe_sw_clr{s,1} + globe_ldw_clr{s,1} - globe_lwup_clr{s,1} ;
end

for s = setdiff(1:18967,del4)
      for i = 1:length(days) % Running for each day

jt = [];
m = [];
idx = [];
G = []; 
if Rin_clr{s,1}(i,1)  - globe_lwtoa_clr{s,1}(i,1) > 3   
jt = (1:(Rin_clr{s,1}(i,1)  - globe_lwtoa_clr{s,1}(i,1)))';

for w = 1:length(jt) % calculate the power and temperature for every possible value of J
    ts(w,1) = ((Rin_clr{s,1}(i,1)  - jt(w,1))/(5.67 * 10^-8)).^0.25;
G(w,1) = (jt(w,1) - du_clr{s,1}(i,1)) * ((ts(w,1)/T_rad_clr{s,1}(i,1)) - 1); % for dissipative heat engine
end
% Numerically solving for dG/dJ = 0

[m,idx] = max(G); % idx represent turbulent flux (J) at max power
mp{s,1}(i,1) = max(G);
turb_clr{s,1}(i,1) = idx;
Tmp_clr{s,1}(i,1) = ((Rin_clr{s,1}(i,1) - turb_clr{s,1}(i,1))/(5.67 * 10^-8)).^0.25;
jt = []; 
m = [];
idx = [];
G = [];
else % Removing the cases that undergoes inversion

idx = 0;
Tmp_clr{s,1}(i,1) = ((Rin_clr{s,1}(i,1) - globe_rnet_clr{s,1}(i,1))/(5.67 * 10^-8)).^0.25;
turb_clr{s,1}(i,1) = globe_rnet_clr{s,1}(i,1);
jt = [];
m = [];
idx = [];
G = [];  
end
    end
    disp(s)
end

% The outputs turb_all and turb_clr represent the optimized turbulent
% fluxes at maximum power for all-sky and clear-sky conditions

% The outputs tmp_all and tmp_clr represent the estimated surface
% temperatures at maximum power for all-sky and clear-sky conditions

% Quantifying cooling and adjusting for clouds
for s = 1:length(grid_mask)
      g_cooling{s,1} = Tmp_clr{s,1} - Tmp_all{s,1};
      g_ta_clear{s,1} = g_tmean_cpc{s,1}(1:5844,1) + g_cooling{s,1};
end


