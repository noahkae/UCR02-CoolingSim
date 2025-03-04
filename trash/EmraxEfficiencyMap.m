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

% create 3D plot
figure;
surf(speedGrid, torqueGrid, efficiencyGrid);

% save to vars
speed = linspace(min(rpm), max(rpm), 100);
torque = linspace(min(torque), max(torque), 100); 
efficiency = efficiencyGrid;

% plot formatting
xlabel('Speed (rpm)');
xlim([0, 5000])
ylabel('Torque (Nm)');
ylim([0,250])
zlabel('Efficiency');
zlim([0,1])
title('EMRAX 228 Efficiency Map');
shading interp;         
colormap winter;