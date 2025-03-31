%F28 Cooling Calculations
clc; clear; 

%% Constants

CpWat = 4.1813; %kj/kg*C
CpStm = 2.08; %kj/kg*C
CpAir = 1.012; %kj/kg*C
InitPress = 101325; %Pa
Rair = 287.058; %J/kg*R
FOS = .66

%% Inputs

TorqueRange = 0:5:200; %Range of torques to calculate cooling load
waterTemps = 25:1:200;
Tamb = 40; %C
FrontalArea = .071/2; %M^2  F26 = .071 **one big radiator double size  F25 = .0565 **one radiator
HeatTransferArea = FrontalArea*46; %m^2 place holder for now reference: F26= 3.3, F25= 2.6
%heat transfer area is calculated to be 46 times the frontal area for the standard fin density
SpeedRatio = .66; % Percent speed of air through radiator relative to car speed  F25 = .52  F26 = .66
HeatTransferCoefficient = .4615; %Must be from empirical testing. .46 was found with f25 and mishimoto radiator. Should be related to construction of radiator only.

CarWaterm = 1; %kg of water in the cooling system

%% Density Change Water, calculation for equation in F27 cooling spreadsheet

% ***not used in f28 stuff, yet***
Twat = 0:1:180; %C
TwatK = Twat+273; %K
RhoWat = -2.86508609805697E-12.*TwatK.^6+6.38960553206707E-09.*TwatK.^5-5.95203372006815E-06.*TwatK.^4+0.00296781391817375.*TwatK.^3-0.838664938806037.*TwatK.^2+127.427758320699.*TwatK-7109.55709155285; %Density of Water Based on Temp [kg/m^3]

%% Density Change Air, calculation for equation in F27 cooling spreadsheet

Tair = 0:1:45; %C
TairK = Tair+273; %K
RhoAir = InitPress./(Rair.*TairK); %Density of Air Based On Temp [kg/m^3]

%% RPM and Speed
gearRatio = 2.65;
WheelCir = pi*18; %in
RevsPerMile = 1/(WheelCir/12/5280);
SpeedArray = 0:1:70; % [MPH]
SpeedArrayMS = SpeedArray.*.44704; % [m/s]
rpmArray = SpeedArray*gearRatio*63360/(WheelCir*60);

% for i = 1:1:size(SpeedArray,2)
%     RPMArray(i)=SpeedArray(i)*RevsPerMile*gearRatio/60;
% end
%% Volumetric Flow of Air

FrontalArea = FrontalArea; %M^2  F26 = .071 **one big radiator double size  F25 = .0565 **one radiator
%changed previously, placeholder
LPSAir = SpeedArrayMS.*FrontalArea*SpeedRatio*1000; %L/s

%% Volumetric Flow of Water

for i = 1:1:size(SpeedArray,2)
    LPSWat(i) = 8/60; % L/s   8L/min to L/s
end

%% Mass Flow of Air

AmbDensity = RhoAir(find(Tair == Tamb)); %kg/m^3
MdotAir = LPSAir./1000.*AmbDensity; %kg/s
for i = 2:1:size(waterTemps,2)
    MdotAir(i,:)=MdotAir(1,:); %kg/s
end

%% Mass Flow of Water
for j = 1:1:size(waterTemps,2)
    WatDensity(j,1) = 997; %kg/m^3
end
MdotWat = LPSWat./1000.*WatDensity; %Mdot Water   kg/s

%% Radiator Cooling 
Ch = MdotWat*CpWat;  %Heat capacity rate water
Cc = MdotAir*CpAir;  %Heat capacity rate air
for i = 1:1:size(SpeedArray,2)
    for j = 1:1:size(waterTemps,2)
        Cmin(j,i) = min(Ch(j,i),Cc(j,i));
        Cmax(j,i) = max(Ch(j,i),Cc(j,i));
    end
end
Cr = Cmin./Cmax;
for i = 1:1:size(waterTemps,2)
    DeltaT(i,1:size(SpeedArray,2)) = waterTemps(i)-Tamb;
end
% for i=1:1:size(SpeedArray,2)
%     DeltaT(:,i) = DeltaT(:,1);
% end

qmax = Cmin.*DeltaT;

NTU = (HeatTransferCoefficient*HeatTransferArea)./Cmin;

for i = 1:1:size(SpeedArray,2)
    for j = 1:1:size(waterTemps,2)
        e(j,i) = 1-exp((-1/Cr(j,i))*(1-exp(-Cr(j,i)*NTU(j,i))));
    end
end
q = qmax.*e;
%Coco's code for radiation
% qRadiation = transpose(.09 * 5.6703*10^-8 * (Th.^4-Tamb^4) * HeatTransferArea);
%
% q = (q + qRadiation).*FOS;
% q(:,1) = qRadiation;
q = q.*FOS;

%% Radiator Cooling - Radiators in parallel, flow/2
Ch = MdotWat*CpWat;  %Heat capacity rate water
Ch2 = Ch/2;
Cc = MdotAir*CpAir;  %Heat capacity rate air
for i = 1:1:size(SpeedArray,2)
    for j = 1:1:size(waterTemps,2)
        Cmin(j,i) = min(Ch2(j,i),Cc(j,i));
        Cmax(j,i) = max(Ch2(j,i),Cc(j,i));
    end
end
Cr = Cmin./Cmax;
for i = 1:1:size(waterTemps,2)
    DeltaT(i,1)=waterTemps(i)-Tamb;
end
for i=1:1:size(SpeedArray,2)
    DeltaT(:,i)=DeltaT(:,1);
end
qmax2 = Cmin.*DeltaT;
NTU = (HeatTransferCoefficient*HeatTransferArea)./Cmin;
for i = 1:1:size(SpeedArray,2)
    for j = 1:1:size(waterTemps,2)
        e(j,i) = 1-exp((-1/Cr(j,i))*(1-exp(-Cr(j,i)*NTU(j,i))));
    end
end
q2 = (qmax2.*e).*FOS;
%% Heat Rejected From Motor and Motor Controller

%*** old, not used right now ***
motorControllerEfficiency = 0.97;

for i = 1:1:size(rpmArray,2)
    for j = 1:1:size(TorqueRange,2)
        Q(j,i) = ((1-calculateEfficiency(rpmArray(i),TorqueRange(j)))...
            *rpmArray(i)*TorqueRange(j)*(motorControllerEfficiency)...
            +(1-motorControllerEfficiency)*rpmArray(i)*TorqueRange(j))/1000/9.5488;
    end
end


% surf(SpeedArray,waterTemps,q)
% xlabel('Speed[MPH]')
% ylabel('Water Temp [C]')
% zlabel('Kilowatts')
% title('Cooling Potential')
% 
% figure
% 
% surf(rpmArray,TorqueRange,Q)
% xlabel('Speed[RPM]')
% ylabel('Torque [NM]')
% zlabel('Kilowatts')
% title('Cooling Load')


%% ************** Endurance Sim **************

%% Import and setup data
%xtime	aps	  speed [km/h]	Torque

data = csvread('HeyFilipUseThisOneModified.csv'); % Import csv file 

%can data files
% time = data(:,1)-data(1,1);
% SimPower = data(:,2)/1000; % kw
% SimRPMs = -data(:,3); %rpm
% SimSpeeds = SimRPMs*pi*18*60/(gearRatio*63360); % mph
% SimTorque = data(:,4); %NM
% actualMotorTemp = data(:,5); %C
% actualBoardTemp = data(:,6); %C

%F27 north autocross [time speed torque]
time = data(:,1)-data(1,1);
SimSpeeds = data(:,3)*0.621371; %mph (originally kpm)
SimTorque = data(:,2)*120; %NM
SimRPMs = SimSpeeds*(gearRatio*63360)/(pi*18*60); %rpm
SimPower = SimTorque.*SimRPMs/9.5488/1000; % kw

%% Running System Temp

OperatingT = zeros(length(SimSpeeds),1);

RunT = Tamb(1,1);

q(:,1) = 0;
q2(:,1) = 0;
timeStep = (time(end) - time(1))/size(time,1);

% SimPower(:,1) = 100;
% SimRPMs(:,1) = 3000;
% SimSpeeds = SimRPMs*pi*18*60/(gearRatio*63360); % mph

for i = 1:1:length(SimSpeeds)
    %Q
    torqueID = abs(round(SimTorque(i)/(TorqueRange(2)-TorqueRange(1))))+1; %Row Coordinate
    SimSpeedID = round(SimSpeeds(i))+1; % Column Coordinate

    overallInefficiency = motorControllerEfficiency*(1-calculateEfficiency(SimRPMs(i),SimTorque(i))) + (1-motorControllerEfficiency);
    SimQadd = SimPower(i)*overallInefficiency/(1-overallInefficiency)+.078; 
    %SimQadd = SimQadd*1.5;
    %SimQadd = 0.11*abs(SimPower(i));
    
    %SimQadd = abs(SimPower(i)*motorControllerEfficiency*.09+SimPower(i)*(1-motorControllerEfficiency)); %constant efficiency
%     SimPower(i)*motorControllerEfficiency*(1-calculateEfficiency(SimRPMs(i),SimTorque(i)))
%     SimPower(i)*(1-motorControllerEfficiency)
    
    %q
    waterTempID = round(RunT + SimQadd/Ch(1,1))-waterTemps(1)+1; %Row Coordinate
    %column coordinate is speed (calculated above)

    %one radiator
    SimQrej = q(waterTempID, SimSpeedID); %Row coordinate

    %%two radiators series
%     SimQrej = q(waterTempID, SimSpeedID);
%     waterTempID = round(RunT + SimQadd/Ch(1,1) - SimQrej/Ch(1,1))-waterTemps(1)+1;%recalculate water temp for the 2nd radiator
%     SimQrej2 = q(waterTempID, SimSpeedID);
%     SimQrej = SimQrej+SimQrej2;
    
    %%two radiators parallel
     %SimQrej = 2*q2(waterTempID, SimSpeedID);

    %%radiator - motor - radiator - motor controller 
%     waterTempID = round(RunT + abs(SimPower(i)*motorControllerEfficiency*(1-calculateEfficiency(SimRPMs(i),SimTorque(i))))/...
%         Ch(1,1))-waterTemps(1)+1; %Water temp ID after motor heat
%     SimQrej = q(waterTempID, SimSpeedID);
%     waterTempID = round(RunT + (SimQadd - SimQrej)/Ch(1,1))-waterTemps(1)+1;%water temp id after motor heat, 1st radiator, motor controller
%     SimQrej2 = q(waterTempID, SimSpeedID);
%     SimQrej = SimQrej+SimQrej2;
    
    
    if i == length(SimSpeeds) %accounts for time end conditions
        SystemDeltaT = (SimQadd-SimQrej)*(time(i)-time(i-1))/(CarWaterm*CpWat);
    else
        SystemDeltaT = (SimQadd-SimQrej)*(time(i+1)-time(i))/(CarWaterm*CpWat);
    end
    
    OperatingT(i,1) = RunT;
    
    RunT = RunT+SystemDeltaT;
    
    QPlot(i) = SimQadd;
    qPlot(i) = SimQrej;
end

% sum = 0;
% for i = 1:1:length(SimSpeeds)
%     sum = sum + qPlot(i)*timeStep;
%     
% end
% sum


% figure
% hold on
% plot(time,QPlot)
% plot(time,qPlot)
% hold off

figure

plot(time,OperatingT)
%legend('Sim Water Temp')
xlabel('Time [s]')
ylabel('Inlet Temp [C]')
ylim([35 55])
%yticks(linspace(0,150,16));

% plot(time,actualMotorTemp,time,actualBoardTemp,time,OperatingT)
% legend('Actual Motor Temp','Actual Board Temp', 'Sim Water Temp')
% xlabel('Time [s]')
% ylabel('Temp (water calculated, driver temp actual) [C]')
% ylim([0 150])
% yticks(linspace(0,150,16));

%Correlation = corr(ActualTmot,OperatingT)^2
