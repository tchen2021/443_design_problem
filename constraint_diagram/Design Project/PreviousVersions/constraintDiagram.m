% Constraint Diagram
clc;clear;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SuperTucano %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumptions 
e= .8; %oswald efficiency factor from range of .7-.85 from wikipedia
h=0;% altitude range [ft]
CD_0=.05; % Random number from agricultural aircraft
%% Input data values
%%%%%% SUPER TOCANO VALUES %%%%%%
S=208.8198621; %ft^2
V_s=121.73333333; % [ft/s]
V_max=469.33333333;
W=11904.962158; % [lb]
P=1600; % [hp] From engine page on wikipedia
AR=6.7; % Aspect Ratio
k=1/(pi*AR*e);
WS=W/S;
S_TOFL=2950; %ft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Standard Atmo Values %%%
[p0, T0, rho0, a0] = atmosphere(0); % Calculate standard atmosphere values at sea level
[p, T, rho, a] = atmosphere(h);% Calculate standard atmosphere values at variable height
% Correct for non standard atmosphere
sigma=rho./rho0;

%% Create Drag Polar
% Calculate Cl_max values at max speed and stall speed
V=V_s:V_max;
% Calc Drag Polar for CD and CL values along velocity
[CD,Cl,Cl_maxL,Cl_max] = dragPolar(W,rho,S,V,V_max,V_s,AR,e,CD_0);

% Plot data
figure(); plot(CD,Cl); hold on; plot(0,CD_0,'r*'); hold on; plot(2*CD_0,sqrt(CD_0/k),'bl*');
grid on; xlabel('CD'); ylabel("CL"); title('Drag Polar - SuperTucano'); 
legend("Drag Polar","CD_0","L/D max",Location="south");

% ASSUMPTIONS - [CD_0, e]

%% Calculate Takeoff Parameters & Graph
%%% ASSUMPTIONS - [.76*CL_MaxL=CL_maxTO]
% The .76% of Cl_MaxL was taken by comparing average mil Patrol
% mean(CL_maxTO)/mean(CL_MaxL). This is just a starting point
CL_maxRatio=.76;
% 1. Calculate Cl_MaxTO from above assumption
Cl_maxTO=CL_maxRatio*Cl_maxL;
% % % % % % % % Cl_maxTO=WS/((T/W)*TOP*S_TOFL*sigma); HOW TO CALC ??????
%%% Calculate TOP (Takeoff Parameters) %%%
TOP=1/( (WS) / ( (P/W) *Cl_maxTO*sigma));

figure(); plot(1/TOP,S_TOFL,'r*'); title("STOFL vs TOP - find best line")
 grid on; xlabel("TOP"); ylabel("STOFL (ft)"); legend("SuperTucano")






% %% Takeoff Distance
% %%% Takeoff Values %%%
% S_TOFL=1500; % [ft] IN ACCORDANCE WITH STOL DEFINITION FOR MILITARY
% % Set W/S range to 100
% WingloadRange=linspace(0,100,101); 
% 
% % Plot the takoff distance
% TW_TO= @(WS) (WS/ (sigma.*Cl_maxTO*TOP*S_TOFL) );
% figure(); plot(WingloadRange,TW_TO(WingloadRange)); hatchedline(WingloadRange,TW_TO(WingloadRange)); hold on;
% % text(73,.43,'Takeoff Distance');  text(65,.15,'Landing Distance');
% % ylim([0,.45]); xlim([0,100]);
% grid on; xlabel("W / S"); ylabel("T / W");



% %% plot the Landing distance
% [MTOWL,ClMaxList] = landingConstraint(S_l);
% plot(MTOWL*ones(1,46),linspace(0,.45,46)); hatchedline(MTOWL*ones(1,46),linspace(0,.45,46),'r-'); hold on;
% 
% %% Graph constraints


