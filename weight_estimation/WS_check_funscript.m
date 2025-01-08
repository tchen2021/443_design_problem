clear
close all
format short
clc

WS = linspace(0,100);
W_TO = [];
for i = 1:length(WS)
    W_TO(i) = Q2_WS_Check(WS(i));
end  
hold on; grid minor
plot(WS,W_TO, 'b', LineWidth=2)
plot(WS,W_TO, 'o', Color='b')
title("W_t_o vs W/S")
xlabel("W/S [lb/ft^2]")
ylabel("W_t_o [lb]")


WTO_vs_WS_turboprop_profile1 = [WS' W_TO'];
%% Weight_fractions_2_turboprop
% rho_lt = 0.001496;       % slugs/ft^3, 15000ft
% Cl_lt = 0.666667;     % CL_Endurance
% Vlt_0 = sqrt( (2*WS) / (rho_lt*Cl_lt) );  % ft/s
% 
% SFC = 0.5; c_j = SFC;
% W = 
% W_pl = 3820;
% Range = [300 300];
% E = 4;
% V_lt = 
% L_Dcr = 15.2805514300000;
% L_Dlt = 13.1768546325000;
% EWF = 0.52;
% RFF = 1;
% initial_guess = 
% p = 1;      % profile
% etap = 0.85;
% % 
% % 
% % 
% [W_TO, W_E, W_F,MFF] = Weight_fractions_2_turboprop(c_j, W, W_pl, Range, E, V_lt, L_Dcr, L_Dlt, etap, EWF, RFF, initial_guess,p)
% % 
% % 

