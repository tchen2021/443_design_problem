% Constraint Diagram - Group 5 - Roman Niemiec
clc;clear;close all;format short;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint Diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pick which loading Method you desire
    % load planeData
    fileLocation="C:\Users\olive\Documents\Aero\Aero433\planeData.xlsx";
    
    planeData = importplaneData(fileLocation, "Sheet1", [2, Inf]);
    
    % Classify Data8uj 
    names=planeData{:,1};
    MTOW=planeData{:,2};
    W_e=planeData{:,3};
    S=planeData{:,4};
    AR=planeData{:,5};
    V_s=planeData{:,6};
    V_max=planeData{:,7};
    P=planeData{:,8};
    S_TOFL=planeData{:,9};
    S_L=planeData{:,10};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumptions  Overview
h=0;% altitude range [ft]
% Assume Cl_MaxTO = Cl_MaxL


%%% Standard Atmo Values %%%
[p0, T0, rho0, a0] = atmosphere(0); % Calculate standard atmosphere values at sea level
[p, T, rho, a] = atmosphere(h);% Calculate standard atmosphere values at variable height
% Correct for non standard atmosphere
sigma=rho./rho0;
hceil = 25000;
[pceil, Tceil, rhoceil, aceil] = atmosphere(hceil);
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
    
    plot(1./TOP(i),S_TOFL(i),'*'); title("STOFL vs TOP - find best line"); hold on;

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




% current objective for cruise constraint

% find a solid AR and oswald efficiency

% justify wing loading numbers

% get CDO from thai sheng chen




% TOverWOA= @(W_S) (C_D_0*q./(W_S)+W_S*1./(q*pi*A*e));


% this wing loading this speed this weight this drag, this sfc Power

% VH = SH* Lh/Sw *C bar =0.4




% %cruise at an altitude
% WshaftverPcruceil = @(W_S) ((rhoceil/rho)^(3/4))*(550*np)*((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A*e)))))*V_cr).^-1;
% WshOverW_Pspeedcruiseceil = WshaftverPcruceil(W_S);

% POverWsuper = (550*np)* ((C_D_0*q./(1.1905e4/208.8199)+(1.1905e4/208.8199)*1./(q*pi*6.7*e))*V_cr).^-1;
% POverWAt6 = (550*np)*((C_D_0*q./(10000/178.68)+(178.68)*1./(q*pi*6.5*e))*V_cr).^-1;  
% POverWOV10= (550*np)*((C_D_0*q./(3.1844e4/290.94)+(3.1844e4/290.94)*1./(q*pi*5.4975*e))*V_cr).^-1;
% 
%  POverWOA1K= (550*np)*((C_D_0*q./(16000/401)+(16000/401)*1./(q*pi*8.7545*e))*V_cr).^-1;
%  POverWKAI= (550*np)*((C_D_0*q./(7308.3/172.3)+(7308.3/172.3)*1./(q*pi*7*e))*V_cr).^-1;
%  POverWArch= (550*np)*((C_D_0*q./(14799.6/401)+(14799.6/401)*1./(q*pi*8.8*e))*V_cr).^-1;

% figure(2)

% Plot take-off constraint
% plot(W_S,WshOverW_Pspeed)
% PshOverW_Speed
% hold on
% 
% plot(1.1905e4/208.8199,POverWsuper,'*')
%  plot(10000/178.68,POverWAt6,'*')
%  plot(3.1844e4/290.94,POverWOV10,'*')
%  plot(3.1844e4/290.94,POverWOA1K,'*')
% plot(7308.3/172.3,POverWKAI,'*')
% plot(14799.6/401,POverWArch,'*')
% hold off

% % Formatting the plot
% grid on;
% xlabel('Wing Loading (W/S) [lb/ft^2]');
% ylabel('Power Loading Ratio (W/P)');
% title('Constraint Diagram');
% legend('Speed Requirement');


% true airspeed and if so at what altitude?
%paul used ceiling altidue

% ceiling requirement
% Given variables
%aircrfats above 25000 feet are suppose to be pressurized

% proppeler at sea level
% for a different altitude

%% speed requirement 

% min cruise speed should be greater than 280 knotts
% values most likely needed q dynamic density at sea level, W_S wing loading, 
% assuming a AR, and e rn
W_S =linspace(0, 100, 60); % lb/ft^2 wing loading
 np = 0.8;% makes most since to use

np1 =0.75;
np2 =0.85; 

V_cr = 472.587;%ft/s

A = 7;%choosing this one
A1 = 6;
A2 = 6.5;


e =0.7; % assumption for all of them

e2 = .75;
e3 = .7;

% e2 =  (1.78*(1-.045*(A)^.68))-0.64;%straight wing formula


q= (1/2)*rho*V_cr^2;%same q for all of them
C_D_0 = 0.020;% from thai sheng chen


%% varied np going with 0.8 average value

PshaftverW = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A*e)))))*V_cr)./(550*np);
PshOverW_speed = PshaftverW(W_S);

PshaftverW_1 = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A*e)))))*V_cr)./(550*np1);
PshOverW_speed_1 = PshaftverW_1(W_S);

PshaftverW_2 = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A*e)))))*V_cr)./(550*np2);
PshOverW_speed_2 = PshaftverW_2(W_S);







%% varying AR goinmg with 7? doesnt vary much at lower wing loading
PshaftverWa1 = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A*e)))))*V_cr)./(550*np);
PshOverW_speeda1 = PshaftverWa1(W_S);

PshaftverW_1a2 = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A1*e)))))*V_cr)./(550*np);
PshOverW_speed_1a2 = PshaftverW_1a2(W_S);

PshaftverW_2a3 = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A2*e)))))*V_cr)./(550*np);
PshOverW_speed_2a3 = PshaftverW_2a3(W_S);


%% varying e goinmg with .7 cant get good oswald like 8 due to cannon wing
PshaftverWa1e = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A*e)))))*V_cr)./(550*np);
PshOverW_speeda1e = PshaftverWa1e(W_S);

PshaftverW_1a2e = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A1*e2)))))*V_cr)./(550*np);
PshOverW_speed_1a2e = PshaftverW_1a2e(W_S);

PshaftverW_2a3e = @(W_S) ((((C_D_0*q./(W_S))+((W_S)*(1./(q*pi*A2*e3)))))*V_cr)./(550*np);
PshOverW_speed_2a3e = PshaftverW_2a3e(W_S);





CL = 0.2; %guess rn need to update
CD = 0.03; %guess rn need to update
sigmaceil = rhoceil/rho; %at sea level
RC = 3250/60; %ft/sec using super tucano need to update
RCP = RC/33000;

%% ceil
% P_wceil= @(W_S) (RCP+ (sqrt(W_S) ./ (19 * (CL^(3/2) ./ CD) * sqrt(sigmaceil))))./(550*np);
% Psh_Wceil = P_wceil(W_S);
% 

%ceiling at 25000 ft
sigmaceil = rhoceil/rho; %at sea level
Vv =  1.667;% ft/s obtiained from gudmunson no clue were they got pg 59
CDmin = .025; %value obtained fromchart in gudmunson pg 59
k = 1./(pi*A*e);
P_wceil = @(W_S) (((Vv./(sqrt(((2/rhoceil).*(W_S)).*sqrt(k/(3*CDmin))))) + (4*sqrt((k*CDmin)/3)))*V_cr)/(550*np);
Psh_Wceil = P_wceil(W_S);






% takoff using ground friction
% ip = 5.75 for sonstant speed prop, 4.60 for fixed pitch props
% Given values (replace these with actual values)
ip = 4.60;%power ,oading
s_TOG = 2275;         % Takeoff ground run distance
k1 =.0376;              % Constant k1
k2 =ip*((1/20)^(1/3)) ;    %proppeller disk loading range from 10-30  
CL_max_TO = 2.0; % Maximum lift coefficient at takeoff
mu_G = .30;     % Ground friction coefficient soft ground worst case

% Calculate (X/W)_TO
P_W_TOground = ( (k1 * (W_S)) / (s_TOG * rho * CL_max_TO) + mu_G + 0.72 * C_D_0 ) / k2;















figure(3)
% Plot take-off constraint
% hold on 
plot(W_S,PshOverW_speed)%speed requi
hold on
plot(W_S,Psh_Wceil)

plot(W_S,P_W_TOground)

hold off

% plot(W_S,WshOverW_Pspeedcruiseceil)
% Formatting the plot
grid on;
xlabel('Wing Loading (W/S) [lb/ft^2]');
ylabel('Power Loading Ratio (P/W)');
title('Constraint Diagram');
legend('speed','ceiling','takeoff');

% h1 = hatchedline(W_S, Wsh_Pceil, 'b');
% h2= hatchedline(W_S, PshOverW_speed, 'r');

% h5 =hatchedline(W_S, WshOverW_Pspeedcruiseceil, 'g');
% plot(1.1905e4/208.8199,POverWsuper,'*')
%  plot(10000/178.68,POverWAt6,'*')
%  plot(3.1844e4/290.94,POverWOV10,'*')
%  plot(3.1844e4/290.94,POverWOA1K,'*')
% plot(7308.3/172.3,POverWKAI,'*')
% plot(14799.6/401,POverWArch,'*')
% xlim([-10 130])
% hold off

figure(4)
% Plot take-off constraint
plot(W_S,Psh_Wceil)
hold on 
plot(W_S,PshOverW_speeda1)%speed requi
plot(W_S,PshOverW_speed_1a2)%speed requi
plot(W_S,PshOverW_speed_2a3)%speed requi

hold off

% plot(W_S,WshOverW_Pspeedcruiseceil)
% Formatting the plot
grid on;
xlabel('Wing Loading (W/S) [lb/ft^2]');
ylabel('Power Loading Ratio (P/W)');
title('Constraint Diagram');
legend('13','6','6.5');




figure(5)
% Plot take-off constraint
% plot(W_S,Wsh_Pceil)
% hold on 
plot(W_S,PshOverW_speeda1e)%speed requi
hold on
plot(W_S,PshOverW_speed_1a2e)%speed requi
plot(W_S,PshOverW_speed_2a3e)%speed requi

hold off

% plot(W_S,WshOverW_Pspeedcruiseceil)
% Formatting the plot
grid on;
xlabel('Wing Loading (W/S) [lb/ft^2]');
ylabel('Power Loading Ratio (P/W)');
title('Constraint Diagram');
legend('.8','.75','.7');


%ceiling a straight line about .05 .1
%cruise speed curving down starting at 0.3 
%FL 250 280kts TAS


















%%jets 

%cruise speed at sea level

T_Wcr= @(W_S) ((C_D_0*q./(W_S)+((W_S)*1./(q*pi*A*e))));
T_Wcruise = T_Wcr(W_S);

%at another altidue T_W= @(W_S) ((rho_0/rhoalt)^(3/4))*((C_D_0*q./(W_S)+(W_S)*1./(q*pi*A*e)));






%ceiling at 25000 ft
Vv =  1.667;% ft/s obtiained from gudmunson no clue were they got pg 59
CDmin = .025; %value obtained fromchart in gudmunson pg 59
k = 1./(pi*A*e);
T_Wce = @(W_S) ((Vv./(sqrt(((2/rhoceil).*(W_S)).*sqrt(k/(3*CDmin))))) + (4*sqrt((k*CDmin)/3)));
T_Wceiling = T_Wce(W_S);



figure(6)
% Plot take-off constraint
plot(W_S, T_Wcruise)
hold on
plot(W_S,T_Wceiling)
 h3 = hatchedline(W_S, T_Wcruise, 'r');
 h4 = hatchedline(W_S, T_Wceiling, 'b');
hold off
% Formatting the plot
grid on;
xlabel('Wing Loading (W/S) [lb/ft^2]');
ylabel('Thrust-Weight Ratio (T/W)');
title('Constraint Diagram');
legend('Ceiling 25000 ft','Speed Requirement at sea level');
% xlim([-10 130])













 %% romans atmo
 function [p, T, rho, a] = atmosphere(h)
% atmosphere(h) function compute p, T, rho and a,  giving the altitude h 
% p : pressure in lb/ft^2
% T : temperature in degree Rankine
% rho : density in slug/ft^3
% a : speed of sound
% h : altitude in ft, h<=104990ft 

%% Constants
Cv = 4290;      % constant volume specific heat for the air
R = 1716;       % ideal gas constant
Cp = R+Cv;      % constant pressure specific heat for air
gamma = Cp/Cv;  % specific heat ratio
T_sl = 518.69;  % sea level temperature
p_sl = 2116.2;  % sea level pressure
%% temperature ratio
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
        
%% pressure density

delta = ( 1-(6.875*(10^-6)).*h ).^5.2561.* ((0 <= h) & (h <= 36089)) +...
        (0.2234*exp(4.806*(10^-5).*(36089-h))).* ((36089 < h )& (h <= 65617)) +...
        (3.174716*(10^-6).*(0.75189+1.0577.*(10^-6).*(h-65617)).^-34.164).* ((65617<h) & (h<=104990));
%% p, T and rho

T = T_sl.*theta;
p = p_sl.*delta;
rho = p./(R.*T);
a = sqrt(gamma*R.*T);

 end




