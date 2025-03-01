% load CSV
data = readtable('LapTime_Power_new__10.xlsx'); 

% define data vectors
time = data{:, 1};
rpm = data{:, 2};
nmtq = data{:, 3};
pwr = data{:, 4};
eta = data{:, 5}/100; % divide by 100 to get percentage /1

% calculate motor heat generated
motor_q = eta.*pwr;

% define controller efficiency
eta_controller = 0.9;
controller_q = eta_controller.*pwr; % calculate controller heat generated
% this is an overestimation, as the controller loss of power should be
% accounted for in motor heat created.

