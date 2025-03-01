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