% % aero 433 assignment 3
% For the same airplane you did the carpet plot assignment. 
% Using a reasonable L/D, EWF and SFC, based on competitive 
% assessment or literature review. Assume that this airplane will be
% a twin-engine airplane.
% 
% Create a constraint diagram for the following requirements:
% 
% Take off distance - filed length -  less than 5000ft at 8000 ft 
% altitude ISA+10C at maximum take-off weight

% Landing Distance - field length - less than 5000ft at sea level 
% ISA+10C at 80% Maximum take-off weight

% Has a climb gradient greater than 0.024 with one engine inoperative 
% (OEI) at clean configuration, no ground effect and V > 1.2Vs

% Maximum speed greater than M0.82 at 33000ft ISA conditions at not
% more than 65% maximum thrust. 
clc
clear

L_D = 14;%range is 13-16
W = 135000;
EWF = 0.45;%chose values from past hw
SFC =.525;
%W/s over T/W
h2 = 8000;
h1 = 0;
CLmax = 2;%from notes
k = -.0000068756; %from gun ch16
Tdelisa = 10;
% Vs = sqrt(((2*W)/S)/rho*CLmax)
% W_S = ((Vs^2)*rho*CLmax)/2

% takeoff requirment?
[p_8000, T_8000, rho_8000, a_8000] = atmospheresendi(h2);
[p_0, T_0, rho_0, a_0] = atmospheresendi(h1);

%standard sea level
T =518.67 * (1 + k * h2);%atmospheric ambient temp from book

rho_8000isa = (1.233 / (T + Tdelisa)) * (1 + k * h2)^5.2561;

sigma = rho_8000isa/rho_0;

W_S =linspace(0, 125, 20); % lb/ft^2;
TOP25 = 120;%from graph in notes
%T/w x= W/s
TOD_TOverW = W_S/(sigma*CLmax*TOP25);


%% landing

SFL = 5000;%ft given in problem
rho_0isa = rho_0/(1+(10)/(T_0)); %slug/ft^3

V_A = sqrt(SFL/0.3);
V_SL = (V_A/1.3)*1.6878;%1.6878 is converting from kts to ft/sec

W_Slanding = (V_SL)^2*rho_0isa*CLmax/2;



% climb gradient 


%climb gradinet > 0.024 with one engine out, clean config, no ground effect, and V > 1.2V_s
A = 9;
e = 0.75;
C_D_0 = 0.0184;

C_LClean = 1.4;
C_DClean = C_D_0+C_LClean^2/(A*pi*e);

LOverDClean  = C_LClean/C_DClean;

CGR = 0.024;
%TW_climb = (N / (N - 1)) * (LD^-1 + CGR)
T_WCGR = 2*(1/(LOverDClean)+CGR);


%% Max Speed 


%must be greater than M=0.82 at 33000 ISA at 65% thrust

[p_33, T_33, rho_33, a_33] = atmospheresendi(33000);
V_fts = 0.82*a_33;
q_82 = (1/2)*rho_33*V_fts^2;
TOverW = @(W_S) (C_D_0*q_82./(W_S)+W_S*1./(q_82*pi*A*e))/0.65;
TOverW_Speed = TOverW(W_S);

figure(1)

% Plot take-off constraint
plot(W_S, TOD_TOverW, 'r-', 'LineWidth', 1.5);
hold on;
% Plot landing constraint as a horizontal line
xline(W_Slanding, 'b--', 'LineWidth', 1.5);
yline(T_WCGR, 'g--','LineWidth', 1.5);
plot(W_S,TOverW_Speed)
hold off
% Formatting the plot
grid on;
xlabel('Wing Loading (W/S) [lb/ft^2]');
ylabel('Thrust-to-Weight Ratio (T/W)');
title('Constraint Diagram');
legend('Take-off Requirement', 'Landing Requirement','Climb Gradient','Speed');
xlim([0, 100])
ylim([0, 0.6])


%% drag polar
C_l = linspace(0,4);

C_D_func = C_D_0+C_l.^2./(A*e*pi);
tanFunc = @(x) 17.1*x;


figure(2)
plot(C_D_func, C_l)
hold on
yline(CLmax)
hold off
%fplot(tanFunc, [0, 0.2], LineWidth=1.5)
xlim([0, 0.4])
ylim([0, 2.5])
xlabel('$C_d$')
ylabel('$C_l$')
grid on


L_Dd = max(C_l./C_D_func)



%% functions
function [p, T, rho, a] = atmospheresendi(h)
% atmosphere(h) function compute p, T, rho and a, giving the altitude h
% p : pressure in lb/ft^2
% T : temperature in degree Rankine
% rho : density in slug/ft^3
% a : speed of sound
% h : altitude in ft, h<=104990ft
%% Constants
Cv = 4290; % constant volume specific heat for the air
R = 1716; % ideal gas constant
Cp = R+Cv; % constant pressure specific heat for air
gamma = Cp/Cv; % specific heat ratio
T_sl = 518.69; % sea level temperature
p_sl = 2116.2; % sea level pressure
%% temperature ratio
if h<=36089
theta = 1-(6.875*10^-6)*h;
elseif h>36089 && h<=65617
theta = 0.75189;
elseif h>65617 && h<=104990
theta = 0.75189+(1.0577*(10^-6))*(h-65617);
else
theta=nan;
warning('altitude exceeded')
end
%% pressure density
delta = ( 1-(6.875*(10^-6)).*h ).^5.2561.* ((0 <= h) & (h <= 36089)) +...
(0.2234*exp(4.806*(10^-5).*(36089-h))).* ((36089 < h )& (h <= 65617)) +...
(3.174716*(10^-6).*(0.75189+1.0577.*(10^-6).*(h-65617)).^-34.164).*((65617<h) & (h<=104990));
%% p, T and rho
T = T_sl.*theta;
p = p_sl.*delta;
rho = p./(R.*T);
a = sqrt(gamma*R.*T);
end



