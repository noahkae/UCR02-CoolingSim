% load CSV
data = readtable('FIltered_Efficiency_Map.csv'); 

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

%--------------------------------------------------------------------------
% load pump curve CSV
pump_curve = readtable('davies-craig-curves.csv'); 

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
radiator_L = 0.24;
radiator_W = 0.012;
radiator_H = 0.12;
tubes_N = 18;
tube_H = 0.00135;
fin_spacing = 0.001;
wall_thickness = 0.0001;
wall_conductivity = 237;
liquid_pipe_D = 0.01;

tube_H_internal = tube_H - 2 * wall_thickness; % [m]
tube_W_internal = radiator_W - 2 * wall_thickness; % [m]

gap_H = (radiator_H - tubes_N * tube_H) / (tubes_N - 1); % [m]
air_area_flow = (tubes_N - 1) * radiator_L * gap_H; % [m^2]

air_area_primary = tubes_N * 2 * (radiator_W + tube_H) * radiator_L; % [m^2]

fins_N = (tubes_N - 1) * radiator_L / fin_spacing;
air_area_fins = 2 * fins_N * radiator_W * gap_H; % [m^2]

thermal_resistance_primary = wall_thickness / air_area_primary / wall_conductivity; % [K/kW]
thermal_resistance_secondary = 1/(55/1000);

liquid_pipe_A = pi * liquid_pipe_D^2 / 4; % [m^2]
