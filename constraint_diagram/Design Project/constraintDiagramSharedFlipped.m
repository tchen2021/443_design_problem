% Constraint Diagram - Group 5 - Roman Niemiec
clc;clear;close all;format short;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint Diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pick which loading Method you desire
    % load planeData
    fileLocation="C:\Users\bravo_4e3\OneDrive - Cal Poly\AERO 435 - Design I\Design Project\planeData.xlsx";
    planeData = importplaneData(fileLocation, "Sheet1", [2, Inf])
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Graph STOFL vs TOP
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
%     TOP(i)=1/(WS / ( (P(i)/W) *Cl_maxTO*sigma));
    TOP(i)=WS*(W/P(i)) /(Cl_maxTO*sigma)
    plot(TOP(i),S_TOFL(i),'*'); title("STOFL vs TOP - find best line"); hold on;

end
 grid on; xlabel("TOP"); ylabel("STOFL (ft)"); legend(names)


%% Plot the Landing distance
for i =1: height(planeData)
    figure (2)
    plot(V_s(i),S_L(i),"*"); title("Landing distance vs Stall speed"); hold on;
end
grid on; xlabel("Stall Speed (ft/s)"); ylabel("S_L (ft)"); legend(names)


%% Create our Constraint Diagram
%Outputs a vector of maximum takeoff wing loadings 
% with corresponding ClMaxes
S_TOFL=1500;
Range=linspace(1,100,100); % range of W/P values
TOP=1/50;
V=400; %ft/s
eff=.82;

%%% CL ITERATIVE
S_L=1500;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]


% STEP 2 - LANDING ROLL RELATIONSHIP

RoskamSpeed2Dist=.265; % Pretty closely verified through our data

V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s
% Cl_maxTO= @(WS) (2*WS)/(rho*knots2ftperS(V_S)^2)
Cl_maxTO=1.8;


%%% Takeoff %%%
% Equation for W/P (Power Loading)
WP= @(WS) ((550*eff)/V) *((WS/ (sigma*Cl_maxTO*TOP*S_TOFL) ));

% %%% Landing %%%
% WeightRatio=.84; % for turboprops
% WS= (rho*Cl_maxL(WS)*V_s(3)^2)/2;

%% speed requirement 
% min cruise speed should be greater than 280 knotts
% values most likely needed q dynamic density at sea level, W_S wing loading, 
% assuming a AR, and e rn

%jets have not done yet
W_S =linspace(0, 100, 50); % lb/ft^2 wing loading

 V_cr = 472.587;%ft/s
 A = 7.5;%assumption from competitor
 e =0.8; % assumption for all of them
 q= (1/2)*rho*V_cr^2;%same q for all of them
C_D_0 = 0.025;%assuming this for all of them can change
TOverW = @(W_S) (C_D_0*q./(W_S)+W_S*1./(q*pi*A*e));%jets
TOverW_Speed = TOverW(W_S);

% proppeler
np = 0.80;%michaels value

WshaftverP = @(W_S) (((C_D_0*q./(W_S)+(W_S)*1./(q*pi*A*e)))*V_cr)./(550*np);
WshOverW_P = WshaftverP(W_S);

% POverWsuper = (550*np)* ((C_D_0*q./(1.1905e4/208.8199)+(1.1905e4/208.8199)*1./(q*pi*6.7*e))*V_cr).^-1;
% POverWAt6 = (550*np)*((C_D_0*q./(10000/178.68)+(178.68)*1./(q*pi*6.5*e))*V_cr).^-1;  
% POverWOV10= (550*np)*((C_D_0*q./(3.1844e4/290.94)+(3.1844e4/290.94)*1./(q*pi*5.4975*e))*V_cr).^-1;
%  
% POverWOA1K= (550*np)*((C_D_0*q./(16000/401)+(16000/401)*1./(q*pi*8.7545*e))*V_cr).^-1;
% POverWKAI= (550*np)*((C_D_0*q./(7308.3/172.3)+(7308.3/172.3)*1./(q*pi*7*e))*V_cr).^-1;
% POverWArch= (550*np)*((C_D_0*q./(14799.6/401)+(14799.6/401)*1./(q*pi*8.8*e))*V_cr).^-1;


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
% p, T and rho
T = T_sl.*theta;
p = p_sl.*delta;
rho = p./(R.*T);
a = sqrt(gamma*R.*T);

%% Ceiling
% Given variables
hceil = 25000;%aircrfats above 25000 feet are suppose to be pressurized
[pceil, Tceil, rhoceil, aceil] = atmosphere(hceil);
CL = 0.2; %guess rn need to update
CD = 0.03; %guess rn need to update
sigmaceil = rhoceil/rho; %at sea level
RC = 3250/60; %ft/min using super tucano need to update
RCP = RC/33000;
W_Pceil= @(W_S) (RCP + (sqrt(W_S) ./ (19 * (CL^(3/2) ./ CD) * sqrt(sigmaceil))))/(np);
Wsh_Pceil = W_Pceil(W_S);



%%% Graph %%%
figure(3)
% Takeoff Distance
plot(Range,WP(Range),'m'); hatchedline(Range,WP(Range),'m',(300-180)*pi/180); hold on;
WeightRatio=.84; % for turboprops

S_L=1500;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]


% STEP 2 - LANDING ROLL RELATIONSHIP

RoskamSpeed2Dist=.265; % Pretty closely verified through our data

V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s


% Landing Distance
Cl_maxL=2.0;
WS1= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio
Cl_maxL=2.4;
WS2= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio
Cl_maxL=2.6;
WS3= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio


plot(WS1*ones(1,46),linspace(0,15,46)); hatchedline(WS1*ones(1,46),linspace(0,15,46),'r-'); hold on;
plot(WS2*ones(1,46),linspace(0,15,46)); hatchedline(WS2*ones(1,46),linspace(0,15,46),'r-'); hold on;
plot(WS3*ones(1,46),linspace(0,15,46)); hatchedline(WS3*ones(1,46),linspace(0,15,46),'r-'); hold on;
% plot(WS*WeightRatio*ones(1,46),linspace(0,15,46)); hatchedline(WS*WeightRatio*ones(1,46),linspace(0,15,46),'r-'); hold on;

S_L=3000;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]


% STEP 2 - LANDING ROLL RELATIONSHIP

RoskamSpeed2Dist=.265; % Pretty closely verified through our data

V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s


% Landing Distance
Cl_maxL=2.0;
WSD1= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio
Cl_maxL=2.4;
WSD2= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio
Cl_maxL=2.6;
WSD3= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio


plot(WSD1*ones(1,46),linspace(0,15,46)); hatchedline(WSD1*ones(1,46),linspace(0,15,46),'b-'); hold on;
plot(WSD2*ones(1,46),linspace(0,15,46)); hatchedline(WSD2*ones(1,46),linspace(0,15,46),'b-'); hold on;
plot(WSD3*ones(1,46),linspace(0,15,46)); hatchedline(WSD3*ones(1,46),linspace(0,15,46),'b-'); hold on;


% Speed
plot(W_S,WshOverW_P); %hatchedline(W_S,WshOverW_P,'b'); hold on;
% Ceiling
plot(W_S,Wsh_Pceil,'g'); hatchedline(W_S,Wsh_Pceil,'g',90/180*pi);hold on

% Graph Parameters
xlim([0,100]); xlabel("W / S"); ylabel("P/ W"); grid on; ylim([0,4]);
% text(30,.43,'Takeoff Distance'); %text(10,.17,'Climb Gradient'); text(2,.4,'Maximum Speed'); text(65,.15,'Landing Distance');
 text(20,15.5,"MTOW"); text(25,14,'WeightRatio')
 %ylim([0,5.0]);
% legend('landing 2.0','landing 2.4','landing 2.6 ','Vcruise','Ceiling')

% %% plot the Landing distance
% [MTOWL,ClMaxList] = landingConstraint(S_l);
% plot(MTOWL*ones(1,46),linspace(0,.45,46)); hatchedline(MTOWL*ones(1,46),linspace(0,.45,46),'r-'); hold on;
% 
% %% Graph constraints


