%--------------------------------------------------------------------------
% constants
rho_air = 1.14; 

%--------------------------------------------------------------------------
% load CSV
data = readtable('data/FIltered_Efficiency_Map.csv'); 

% extract cols
rpm = data{:, 1};     
torque = data{:, 2};    
efficiency = data{:, 3};

% combine rpm and torque into a matrix 
inputData = [rpm, torque];

% fit 5th-degree polynomial
fitModel = polyfitn(inputData, efficiency, 5);

% generate grid for rpm and torque
[speedGrid, torqueGrid] = meshgrid(linspace(min(rpm), max(rpm), 100), ...
                                   linspace(min(torque), max(torque), 100));

% predict efficiency
efficiencyGrid = polyvaln(fitModel, [speedGrid(:), torqueGrid(:)]);
efficiencyGrid = reshape(efficiencyGrid, size(speedGrid));

% set max efficiency to 96%
efficiencyGrid = min(efficiencyGrid, 0.96);

% save to vars
speed = linspace(min(rpm), max(rpm), 100);
torque = linspace(min(torque), max(torque), 100); 
efficiency = efficiencyGrid;

figure;
surf(speedGrid, torqueGrid, efficiencyGrid, 'EdgeColor', 'none');
colormap('winter'); % MATLAB's blue-based colormap
colorbar;
xlabel('RPM');
ylabel('Torque (Nm)');
zlabel('Efficiency');
zlim([0,1]);
title('Efficiency Map');
view(135, 30); % Adjust the viewing angle

%--------------------------------------------------------------------------
% load pump curve CSV
pump_curve = readtable('data/davies-craig-curves.csv'); 

% save EBP40 24V data
ebp40_24v_lpm = fit_and_interpolate(pump_curve{:, 1});
ebp40_24v_head = fit_and_interpolate(pump_curve{:, 2});

% save EBP40 12V data - shit is wack
%ebp40_12v_lpm = fit_and_interpolate(pump_curve{:, 3});
%ebp40_12v_head = fit_and_interpolate(pump_curve{:, 4});

% save EBP25 data
ebp25_lpm = fit_and_interpolate(pump_curve{:, 5});
ebp25_head = fit_and_interpolate(pump_curve{:, 6});

% save EBP23 data
ebp23_lpm = fit_and_interpolate(pump_curve{:, 7});
ebp23_head = fit_and_interpolate(pump_curve{:, 8});

%--------------------------------------------------------------------------
% plot head vs. LPM
figure; hold on; grid on;
plot(ebp40_24v_lpm, ebp40_24v_head, '-o', 'DisplayName', 'EBP40 24V');
%plot(ebp40_12v_lpm, ebp40_12v_head, '-o', 'DisplayName', 'EBP40 12V');
plot(ebp25_lpm, ebp25_head, '-o', 'DisplayName', 'EBP25');
plot(ebp23_lpm, ebp23_head, '-o', 'DisplayName', 'EBP23');

xlabel('Flow Rate (LPM)');
ylabel('Head (m)');
title('Pump Head vs. Flow Rate');
legend;

%--------------------------------------------------------------------------
function y_out = fit_and_interpolate(x)
    % remove NaNs
    x_clean = x(~isnan(x));
    
    % generate an index vector
    x_idx = 1:length(x_clean);
    
    % fit a polynomial
    p = polyfit(x_idx, x_clean, 10);
    
    % interpolate to evenly spaced points
    x_interp = linspace(1, length(x_clean), 50);
    y_out = polyval(p, x_interp);
end

%--------------------------------------------------------------------------
% Aliexpress rad dimensions and properties
ali_rad_L = 0.24;
ali_rad_W = 0.017;
ali_rad_H = 0.120;
ali_rad_N = 18;
ali_rad_tube_H = 0.0013;
ali_rad_fin_spacing = 0.00075;
ali_rad_wall_thickness = 0.0001;
ali_rad_conductivity = 237;

% function call of Ali rad
ali_results = generate_rad_dims(ali_rad_H, ali_rad_W, ali_rad_L, ali_rad_N, ali_rad_tube_H, ali_rad_wall_thickness, ali_rad_fin_spacing, ali_rad_conductivity, 0.016);
ali_tube_H_internal = ali_results.tube_H_internal;
ali_tube_W_internal = ali_results.tube_W_internal;
ali_gap_H = ali_results.gap_H;
ali_air_area_flow = ali_results.air_area_flow;
ali_air_area_primary = ali_results.air_area_primary;
ali_fins_N = ali_results.fins_N;
ali_air_area_fins = ali_results.air_area_fins;
ali_thermal_resistance_primary = ali_results.thermal_resistance_primary;

% broken mishimoto rad dimensions
mishi_rad_L = 0.24;
mishi_rad_W = 0.04;
mishi_rad_H = 0.127;
mishi_rad_tubes_N = 12;
mishi_rad_tube_H = 0.0016;
mishi_rad_fin_spacing = 0.0013;
mishi_rad_wall_thickness = 0.0003;
mishi_rad_conductivity = 237;

% function call of mishi rad
mishi_results = generate_rad_dims(mishi_rad_H, mishi_rad_W, mishi_rad_L, mishi_rad_tubes_N, mishi_rad_tube_H, mishi_rad_wall_thickness, mishi_rad_fin_spacing, mishi_rad_conductivity, 0.016);
mishi_tube_H_internal = mishi_results.tube_H_internal;
mishi_tube_W_internal = mishi_results.tube_W_internal;
mishi_gap_H = mishi_results.gap_H;
mishi_air_area_flow = mishi_results.air_area_flow;
mishi_air_area_primary = mishi_results.air_area_primary;
mishi_fins_N = mishi_results.fins_N;
mishi_air_area_fins = mishi_results.air_area_fins;
mishi_thermal_resistance_primary = mishi_results.thermal_resistance_primary;

% Current rad dimensions
current_rad_L = 0.20;
current_rad_W = 0.04;
current_rad_H = 0.20;
current_rad_tubes_N = 19;
current_rad_tube_H = 0.002;
current_rad_fin_spacing = 0.001;
current_rad_wall_thickness = 0.0003;
current_rad_conductivity = 237;

% function call of current rad
current_results = generate_rad_dims(current_rad_H, current_rad_W, current_rad_L, current_rad_tubes_N, current_rad_tube_H, current_rad_wall_thickness, current_rad_fin_spacing, current_rad_conductivity, 0.016);
current_tube_H_internal = current_results.tube_H_internal;
current_tube_W_internal = current_results.tube_W_internal;
current_gap_H = current_results.gap_H;
current_air_area_flow = current_results.air_area_flow;
current_air_area_primary = current_results.air_area_primary;
current_fins_N = current_results.fins_N;
current_air_area_fins = current_results.air_area_fins;
current_thermal_resistance_primary = current_results.thermal_resistance_primary;

% mishi YXR700 Rhino rad guess values
yxr_rad_L = 0.378968;
yxr_rad_W = 0.04;
yxr_rad_H = 0.258064;
yxr_rad_tubes_N = 23;
yxr_rad_tube_H = 0.0014;
yxr_rad_fin_spacing = 0.001;
yxr_rad_wall_thickness = 0.0003;
yxr_rad_conductivity = 237;
yxr_pipe_D = 0.016; %not real

% function call of yxr rad
yxr_results = generate_rad_dims(yxr_rad_H, yxr_rad_W, yxr_rad_L, yxr_rad_tubes_N, yxr_rad_tube_H, yxr_rad_wall_thickness, yxr_rad_fin_spacing, yxr_rad_conductivity, yxr_pipe_D);
yxr_tube_H_internal = yxr_results.tube_H_internal;
yxr_tube_W_internal = yxr_results.tube_W_internal;
yxr_gap_H = yxr_results.gap_H;
yxr_air_area_flow = yxr_results.air_area_flow;
yxr_air_area_primary = yxr_results.air_area_primary;
yxr_fins_N = yxr_results.fins_N;
yxr_air_area_fins = yxr_results.air_area_fins;
yxr_thermal_resistance_primary = yxr_results.thermal_resistance_primary;
yxr_pipe_A = yxr_results.liquid_pipe_A;

% mishi crf450x rad guess values
crf_rad_L = 0.23;
crf_rad_W = 0.034;
crf_rad_H = 0.128;
crf_rad_tubes_N = 11;
crf_rad_tube_H = 0.0014;
crf_rad_fin_spacing = 0.001;
crf_rad_wall_thickness = 0.0003;
crf_rad_conductivity = 237;
crf_pipe_D = 0.016;

% function call of crf rad
crf_results = generate_rad_dims(crf_rad_H, crf_rad_W, crf_rad_L, crf_rad_tubes_N, crf_rad_tube_H, crf_rad_wall_thickness, crf_rad_fin_spacing, crf_rad_conductivity, crf_pipe_D);
crf_tube_H_internal = crf_results.tube_H_internal;
crf_tube_W_internal = crf_results.tube_W_internal;
crf_gap_H = crf_results.gap_H;
crf_air_area_flow = crf_results.air_area_flow;
crf_air_area_primary = crf_results.air_area_primary;
crf_fins_N = crf_results.fins_N;
crf_air_area_fins = crf_results.air_area_fins;
crf_thermal_resistance_primary = crf_results.thermal_resistance_primary;
crf_pipe_A = crf_results.liquid_pipe_A;

% Current Rad area for airflow
current_rad_area = crf_air_area_flow;

%--------------------------------------------------------------------------
% function to generate radiator simulink representation variables
function results = generate_rad_dims(radiator_H, radiator_W, radiator_L, tubes_N, tube_H, wall_thickness, fin_spacing, wall_conductivity, liquid_pipe_D)
    % internal tube dimensions
    tube_H_internal = tube_H - 2 * wall_thickness; % [m]
    tube_W_internal = radiator_W - 2 * wall_thickness; % [m]

    % air gap height
    gap_H = (radiator_H - tubes_N * tube_H) / (tubes_N - 1); % [m]
    
    % air flow area
    air_area_flow = (tubes_N - 1) * radiator_L * gap_H; % [m^2]

    % primary heat transfer area
    air_area_primary = tubes_N * 2 * (radiator_W + tube_H) * radiator_L; % [m^2]

    % number of fins
    fins_N = (tubes_N - 1) * radiator_L / fin_spacing;
    
    % finned area
    air_area_fins = 2 * fins_N * radiator_W * gap_H; % [m^2]

    % thermal resistance
    thermal_resistance_primary = wall_thickness / (air_area_primary * wall_conductivity); % [K/kW]

    % pipe area
    liquid_pipe_A = pi * liquid_pipe_D^2 / 4; % [m^2]

    
    % store results in a struct
    results = struct('tube_H_internal', tube_H_internal, ...
                     'tube_W_internal', tube_W_internal, ...
                     'gap_H', gap_H, ...
                     'air_area_flow', air_area_flow, ...
                     'air_area_primary', air_area_primary, ...
                     'fins_N', fins_N, ...
                     'air_area_fins', air_area_fins, ...
                     'thermal_resistance_primary', thermal_resistance_primary, ...
                     'liquid_pipe_A', liquid_pipe_A);
end

%--------------------------------------------------------------------------
% load fan params
data = readtable('data/13.8VFan.csv');
data{:, :} = max(data{:, :}, 0.01);
fan_flow_rate = data{:, 1}./60;     
fan_pressure = data{:, 2}; 
fan_rpm = 6400;
fan_diam = 0.12;
fan_area = pi * fan_diam^2 / 4;
figure;
plot(fan_flow_rate, fan_pressure, '-o', 'DisplayName', '13.8V fan');

xlabel('Flow Rate [m3/s]');
ylabel('Pressure Drop [Pa]');
title('Fan Curve');
legend;
