clc; clearvars; close all

LOverD = 14;
SFC = 0.525; %lb/lb/hr
EWF = 0.45;

C_L_max = 2.2;

WOverS = 0:1:150;

%we are making a constraint diagram of W/S vs T/W (weight, wing area, thrust)


%for take off distance use the take off curve fit equation

% Create a constraint diagram for the following requirements:
% 
%     Take off distance - field length -  less than 5000ft at 8000 ft altitude ISA+10C at maximum take-off weight
%     Landing Distance - field length - less than 5000ft at sea level ISA+10C at 80% Maximum take-off weight
%     Has a climb gradient greater than 0.024 with one engine inoperative (OEI) at clean configuration, no ground effect and V > 1.2Vs
%     Maximum speed greater than M0.82 at 33000ft ISA conditions at not more than 65% maximum thrust. 

% %% Take Off Distance

% S_TOFL \propto (W/S)/(sigma*C_L_max_TO*T/W))

[p_0, T_0, rho_0, a_0] = atmosphereImperial(0);
[p_8000, T_8000, rho_8000, a_8000] = atmosphereImperial(8000);

rho_8000_acc = rho_8000/(1+(10)/(T_8000+273.15));

sigma = rho_8000_acc/rho_0;

TOP25 = 133.3; %based off of roskam figure 3.7

TOD_TOverW = WOverS/(sigma*C_L_max*TOP25);

figure()
hold on
p(1) = plot([1000; 1100], [1000; 1100], 'r');
p(2) = plot([1000; 1100], [1000; 1100], 'b');
p(3) = plot([1000; 1100], [1000; 1100], 'g');
p(4) = plot([1000; 1100], [1000; 1100], 'm');
hatchedline(WOverS, TOD_TOverW, 'r-', 90*pi/180, 1, -0.014, 0.019)
grid on
xlabel('W/S, $\mathrm{\left[\frac{lbf}{ft^2}\right]}$')
ylabel('T/W')
% text(50, 0.7, 'Take off distance')
% text(50, 0.65, '5000 ft feild')
% text(50, 0.6, '8000 ft altitude ISA+10C')
fontsize(24, 'points')
xlim([0, 100])
ylim([0, 0.6])

% %% Landing Distance

W_LOverW_TO = 0.8;
S_FL = 5000;

%I NEED RHO

[~, T_0, rho_0, ~] = atmosphereImperial(0);

rho_0_nISA = rho_0/(1+(18)/(T_0)); %slug/ft^3

V_A = sqrt(S_FL/0.3);
V_SL = V_A/1.3;

WOverS_landing = ((V_SL*1.6878)^2)*rho_0_nISA*C_L_max/2/W_LOverW_TO;

hatchedline(ones([2,1])*WOverS_landing, [0, 1], 'b-')%, 60*pi/180, 1, -0.03, 0.02)

% %% Climb Gradient

%climb gradinet > 0.024 with one engine out, clean config, no ground effect, and V > 1.2V_s

% C_LClean = 1.4;
% C_DClean = 0.0184+C_LClean^2/26.7;
% 
% LOverDClean  = C_LClean/C_DClean;

TOverW_CGR = 2*(1/(LOverD)+0.024);

hatchedline([0, 100], ones([2,1])*TOverW_CGR, 'g-')%, 60*pi/180, 1, -0.03, 0.02)

% Max Speed 

A = 9;
e = 0.75;
C_D_0 = 0.0184;

%must be greater than M=0.82 at 33000 ISA at 65% thrust

[p_33, T_33, rho_33, a_33] = atmosphereImperial(33000);
V_fts = 0.82*a_33;
q_82 = 1/2*rho_33*V_fts^2;

TOverW = @(WOverS) (C_D_0*q_82./(WOverS)+WOverS*1./(q_82*pi*A*e))/0.65;

%C_D_0, A, e

TOverW_Speed = TOverW(WOverS(2:end));

hatchedline(WOverS(2:end), TOverW_Speed, 'm-')

legend(p(1:4), 'Take Off', 'Landing', 'Climb Gradient', 'Max Speed', Location='best')

%% drag polar

C_D_func = @(C_l) C_D_0+C_l.^2./(A*e*pi)
tanFunc = @(x) 17.1*x

space = linspace(0,3)

figure()
plot(C_D_func(space), space)
hold on
yline(C_L_max)
fplot(tanFunc, [0, 0.2], LineWidth=1.5)
xlim([0, 0.4])
ylim([0, 2.5])
xlabel('$C_d$')
ylabel('$C_l$')
grid on
fontsize(20, 'points')

%% Functions

function [p, T, rho, a] = atmosphereImperial(h)
% atmosphere(h) function compute p, T, rho and a,  giving the altitude h 
% p : pressure in lb/ft^2
% T : temperature in degree Rankine
% rho : density in slug/ft^3
% a : speed of sound in ft/s
% h : altitude in ft, h<=104990ft 

    % %% Constants
    Cv = 4290;      % constant volume specific heat for the air
    R = 1716;       % ideal gas constant
    Cp = R+Cv;      % constant pressure specific heat for air
    gamma = Cp/Cv;  % specific heat ratio
    T_sl = 518.69;  % sea level temperature
    p_sl = 2116.2;  % sea level pressure
    % %% temperature ratio
    if  h<=36089
        theta = 1-(6.875*10^-6)*h;
    elseif h>36089 && h<=65617
        theta = 0.75189;
    elseif h>65617 && h<=104990
        theta = 0.75189+(1.0577*(10^-6))*(h-65617);
    else
        theta=nan;
        warning('altitude exceeded')
    end
            
    % %% pressure density
    
    delta = ( 1-(6.875*(10^-6)).*h ).^5.2561.* ((0 <= h) & (h <= 36089)) +...
            (0.2234*exp(4.806*(10^-5).*(36089-h))).* ((36089 < h )& (h <= 65617)) +...
            (3.174716*(10^-6).*(0.75189+1.0577.*(10^-6).*(h-65617)).^-34.164).* ((65617<h) & (h<=104990));
    % %% p, T and rho
    
    T = T_sl.*theta;
    p = p_sl.*delta;
    rho = p./(R.*T);
    a = sqrt(gamma*R.*T);
end




%% Functions

function [p, T, rho, a] = atmosphereImperial(h)
% atmosphere(h) function compute p, T, rho and a,  giving the altitude h 
% p : pressure in lb/ft^2
% T : temperature in degree Rankine
% rho : density in slug/ft^3
% a : speed of sound in ft/s
% h : altitude in ft, h<=104990ft 

    % %% Constants
    Cv = 4290;      % constant volume specific heat for the air
    R = 1716;       % ideal gas constant
    Cp = R+Cv;      % constant pressure specific heat for air
    gamma = Cp/Cv;  % specific heat ratio
    T_sl = 518.69;  % sea level temperature
    p_sl = 2116.2;  % sea level pressure
    % %% temperature ratio
    if  h<=36089
        theta = 1-(6.875*10^-6)*h;
    elseif h>36089 && h<=65617
        theta = 0.75189;
    elseif h>65617 && h<=104990
        theta = 0.75189+(1.0577*(10^-6))*(h-65617);
    else
        theta=nan;
        warning('altitude exceeded')
    end
            
    % %% pressure density
    
    delta = ( 1-(6.875*(10^-6)).*h ).^5.2561.* ((0 <= h) & (h <= 36089)) +...
            (0.2234*exp(4.806*(10^-5).*(36089-h))).* ((36089 < h )& (h <= 65617)) +...
            (3.174716*(10^-6).*(0.75189+1.0577.*(10^-6).*(h-65617)).^-34.164).* ((65617<h) & (h<=104990));
    % %% p, T and rho
    
    T = T_sl.*theta;
    p = p_sl.*delta;
    rho = p./(R.*T);
    a = sqrt(gamma*R.*T);
end


