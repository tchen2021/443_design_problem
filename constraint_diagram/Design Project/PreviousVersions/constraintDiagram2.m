% Constraint Diagram - Group 5 - Roman Niemiec
clc;clear;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint Diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%
load planeData
    names=planeData{:,1};
    MTOW=planeData{:,2};
    S=planeData{:,3};
    AR=planeData{:,4};
    V_s=planeData{:,5};
    V_max=planeData{:,6};
    P=planeData{:,7};
    S_TOFL=planeData{:,8};
    S_L=planeData{:,9};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumptions  Overview
h=0;% altitude range [ft]
% Assume Cl_MaxTO = Cl_MaxL


%%% Standard Atmo Values %%%
[p0, T0, rho0, a0] = atmosphere(0); % Calculate standard atmosphere values at sea level
[p, T, rho, a] = atmosphere(h);% Calculate standard atmosphere values at variable height
% Correct for non standard atmosphere
sigma=rho./rho0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i =1: height(planeData)
    % Quick Conversions
    W=MTOW(i);
    WS=W/S(i);

    %% Calculate Takeoff Parameters & Graph - Got rid of assumed ratio
    %%% Assumption - Assume Cl_MaxTO = Cl_MaxL
    % Calculate Cl_maxL from stall speed
    Cl_maxL=(2*W)/(rho*S(i)*V_s(i)^2);
    % CL_max TO 
    Cl_maxTO=Cl_maxL;
    %%% Calculate TOP (Takeoff Parameters) %%%
    TOP(i)=1/(WS / ( (P(i)/W) *Cl_maxTO*sigma));
    
    plot(TOP(i),S_TOFL(i),'*'); title("STOFL vs TOP - find best line"); hold on;

end

 grid on; xlabel("TOP"); ylabel("STOFL (ft)"); legend(names)


%% Plot the Landing distancea
% plot(S_L,V_s.^2); title("Landing distance vs Stall speed squared");
% grid on; xlabel("Stall Speed (ft/s)"); ylabel("S_L (ft)"); legend(names)
% 


%Outputs a vector of maximum takeoff wing loadings 
% with corresponding ClMaxes





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


