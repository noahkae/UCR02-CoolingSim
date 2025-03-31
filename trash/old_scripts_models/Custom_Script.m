%% Custom UCR-02 Cooling Script

%% Loading Emrax efficiency map into a polynomial
% load CSV
data = readtable('data/FIltered_Efficiency_Map.csv'); 

% extract cols
eff_rpm = data{:, 1};     
eff_torque = data{:, 2};    
efficiency = data{:, 3};

% combine rpm and torque into a matrix 
inputData = [eff_rpm, eff_torque];

% fit 5th-degree polynomial
fitModel = polyfitn(inputData, efficiency, 5);

% generate grid for rpm and torque
[speedGrid, torqueGrid] = meshgrid(linspace(min(eff_rpm), max(eff_rpm), 100), ...
                                   linspace(min(eff_torque), max(eff_torque), 100));

% predict efficiency
efficiencyGrid = polyvaln(fitModel, [speedGrid(:), torqueGrid(:)]);
efficiencyGrid = reshape(efficiencyGrid, size(speedGrid));

% set max efficiency to 96%
efficiencyGrid = min(efficiencyGrid, 0.96);

% save to vars
speed = linspace(min(eff_rpm), max(eff_rpm), 100);
eff_torque = linspace(min(eff_torque), max(eff_torque), 100); 
efficiency = efficiencyGrid;

% plot
figure;
surf(speedGrid, torqueGrid, efficiencyGrid, 'EdgeColor', 'none');
colormap('winter');
colorbar;
xlabel('RPM');
ylabel('Torque (Nm)');
zlabel('Efficiency');
zlim([0,1]);
title('Efficiency Map');
view(135, 30); 

%% Load laptime data from sim
% load csv
lapdata = readtable('data/extra_shortened_velocity.csv');

% load cols into variables
lap_time = lapdata{:, 1};   %[s]
lap_rpm = lapdata{:, 2};    %[rpm]
lap_tq = lapdata{:, 3};     %[N-m]
lap_pwr = lapdata{:, 4};    %[W]
lap_speed = lapdata{:, 5};  %[m/s]

%% Constant Values and environment variables
C_p_h2o = 4.18;     %[kJ/kgK]
C_p_air = 1.006;    %[kJ/kgK]
P_atm = 101.325;    %[kPa]
T_atm = 35;         %[°C]

%% Radiator Parameters
% Aliexpress rad dimensions and properties
ali_rad_L = 0.24;                   %[m]
ali_rad_W = 0.017;                  %[m]
ali_rad_H = 0.120;                  %[m]
ali_rad_N = 18;                     %[Dimensionless]
ali_rad_tube_H = 0.0013;            %[m]
ali_rad_fin_spacing = 0.00075;      %[m]
ali_rad_wall_thickness = 0.0001;    %[m]
ali_rad_conductivity = 237;         %[W/mK]
air_area_primary = ali_rad_N * 2 * (ali_rad_W + ali_rad_H) * ali_rad_L; % [m^2]
gap_H = (ali_rad_H - ali_rad_N * ali_rad_tube_H) / (ali_rad_N - 1); % [m]
air_area_fins = 2 * ali_rad_N * ali_rad_W * gap_H; % [m^2]
air_area_total = air_area_fins + air_area_primary; %[m^2]
h_coeff = 55/air_area_total;