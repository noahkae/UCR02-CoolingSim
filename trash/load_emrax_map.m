data = readtable('FIltered_Efficiency_Map.csv'); 
% extract cols
rpm = data{:, 1};     
torque = data{:, 2};    
efficiency = data{:, 3};