clc; clear; close all;
%% Flight angle of attack calculator

% defining inputs from drag buildup
K_W = 0.21003;
K_B = 1.1306;
Se_S = 0.8;
ST_S = 0.8;
alpha_0W = deg2rad(-7);
alpha_0B = 0;
D = 6.1603;
K = 0.73355;
C_DC = 1.1803;
eta = 0.58894;
S_P_x0 = 21.9360;
C_L_alpha = 4.7598;
S = 224; 
i_w = deg2rad(2);
%C_LWB = 0.4283;
Vbar = 

%% calculating C_LWB
W = [];
V = [];
C_L = 2*W/(S*rho*V^2);
C_LT = (1/Vbar) * (C_M0 + C_L * (h-h0));
C_LWB = C_L - (ST_S) * C_LT;

%% Calculating the aoa
syms alpha

d2 = (K_W-K_B) * C_L_alpha * alpha_0W * Se_S;
d1 = (K_W-K_B) * C_L_alpha * alpha * Se_S;
c3 = -(alpha_0B/S) * ((K*pi*D^2/2) + eta*C_DC * (-alpha_0B) * S_P_x0);
c2 = (alpha/S) * ((K*pi*D^2/2) - 2*eta*C_DC * (alpha_0B) * S_P_x0);
e2 = c3 + d2 + d1 + i_w;
e1 = c2 + d1;
c1 = (1/S) * (eta*C_DC*S_P_x0);
eqn = C_LWB == c1 * alpha^2 + e1 * alpha + e2;

Sol = solve(eqn, alpha)
disp(rad2deg(Sol))