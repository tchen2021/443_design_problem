% Constraint Diagram
clc;clear;close all;


% Assume Roskan pg 91 transport jet Cl_max
Cl_max= 1.5; % 1.2-1.8
Cl_maxTO=1.9; % range of 1.6-2.2
% Cl_maxTO=1.6:.2:2.2; % CL MAX TO FROM 1.6-2.2 in iterations of .2
Cl_maxL=2.6; % range of 1.8-2.8


%% Input airplane values (OV-10D)

% data=[]
% MTOW= 1;



%%% Given Values from problem statement %%%
% Ranged Values
WingloadRange=linspace(0,100,101); % Set W/S range to 100
h=0; % altitude range [ft]

%%% Takeoff Values %%%
S_TOFL=2800; % [ft] at MTOW
TOP25=1/37.5; % sourced from Roskam pg 




% Calculate standard atmosphere values at sea level
[p0, T0, rho0, a0] = atmosphere(0);
% Calculate standard atmosphere values at variable height
[p, T, rho, a] = atmosphere(h);

% Correct for non standard atmosphere
sigma=rho./rho0;



%% Plot the takoff distance
TW_TO= @(WS) (WS/ (sigma.*Cl_maxTO*TOP25*S_TOFL) );
plot(WingloadRange,TW_TO(WingloadRange)); hatchedline(WingloadRange,TW_TO(WingloadRange)); hold on;

%% plot the Landing distance
[MTOWL,ClMaxList] = landingConstraint(S_l);
plot(MTOWL*ones(1,46),linspace(0,.45,46)); hatchedline(MTOWL*ones(1,46),linspace(0,.45,46),'r-'); hold on;

%% Graph constraints
grid on; xlabel("W / S"); ylabel("T / W");
text(73,.43,'Takeoff Distance');  text(65,.15,'Landing Distance');
ylim([0,.45]); xlim([0,100]);
