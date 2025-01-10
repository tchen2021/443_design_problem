% Constraint Diagram - Group 5 - Roman Niemiec
clc;clear;close all;format short;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint Diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pick which loading Method you desire
    % load planeData
    fileLocation="C:\Users\bravo_4e3\OneDrive - Cal Poly\AERO 435 - Design I\Design Project\planeData.xlsx";
    planeData = importplaneData(fileLocation, "Sheet1", [2, Inf]);
    jetplaneData = importplaneData(fileLocation, "Jet", [2, Inf]);
    %%% Classify Data  %%%
    %Turboprop
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
    %turbojet
    namesJ=jetplaneData{:,1};
    MTOWJ=jetplaneData{:,2};
    W_eJ=jetplaneData{:,3};
    SJ=jetplaneData{:,4};
    ARJ=jetplaneData{:,5};
    V_sJ=jetplaneData{:,6};
    V_maxJ=jetplaneData{:,7};
    Thrust=jetplaneData{:,8};
    S_TOFLJ=jetplaneData{:,9};
    S_LJ=jetplaneData{:,10};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumptions  Overview
h=0;% altitude range [ft]

%%% Standard Atmo Values %%%
[p0, T0, rho0, a0] = atmosphere(0); % Calculate standard atmosphere values at sea level
[p, T, rho, a] = atmosphere(h);% Calculate standard atmosphere values at variable height
% Correct for non standard atmosphere
sigma=rho./rho0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Takeoff %%%
S_TOFL=2275;
Range=linspace(1,100,100); % range of W/P values
V=400; %ft/s
eff=.82;
Cl_maxTO=1.8; %%% ASSUMPTION
TOPchange= @(TOPtest) 8.134*TOPtest + .0149*TOPtest^2 - S_TOFL;
TOP=fzero(TOPchange,.1);

% Equation for W/P (Power Loading)
WP= @(WS) WS/(sigma*Cl_maxTO*TOP);
Cl_maxTO=2.0; %%% ASSUMPTION

S_TOFL=3000;
TOPchange= @(TOPtest) 8.134*TOPtest + .0149*TOPtest^2 - S_TOFL;
TOP=fzero(TOPchange,.1);
WP1= @(WS) WS/(sigma*Cl_maxTO*TOP);

%% landing Calcs
WeightRatio=.84;
%1500 ft case
% Landing Distance
S_L=2320;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]
% STEP 2 - LANDING ROLL RELATIONSHIP
RoskamSpeed2Dist=.265; % Pretty closely verified through our data
V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s
Cl_maxL=2.4;
WS= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio;


%% speed requirement 
% min cruise speed should be greater than 280 knotts
% values most likely needed q dynamic density at sea level, W_S wing loading, 
% assuming a AR, and e rn

%jets have not done yet
W_S =linspace(1, 100, 100); % lb/ft^2 wing loading

 V_cr = 472.587;%ft/s
 A = 6.6;%assumption from competitor - average of comp assesment
 e =0.8; % assumption for all of them
 q= (1/2)*rho*V_cr^2;%same q for all of them
C_D_0 = 0.025;%assuming this for all of them can change
TOverW = @(W_S) (C_D_0*q./(W_S)+W_S*1./(q*pi*A*e));%jets
TOverW_Speed = TOverW(W_S);

% proppeler
np = 0.80;%michaels value

WshaftverP = @(W_S) (((C_D_0*q./(W_S)+(W_S)*1./(q*pi*A*e)))*V_cr)./(550*np);
WshOverW_P = WshaftverP(W_S);


%% Ceiling
% Given variables
hceil = 25000;%aircrfats above 25000 feet are suppose to be pressurized
[pceil, Tceil, rhoceil, aceil] = atmosphere(hceil);
%ceiling at 25000 ft
sigmaceil = rhoceil/rho; %at sea level
Vv =  1.667;% ft/s obtiained from gudmunson no clue were they got pg 59
CDmin = .025; %value obtained fromchart in gudmunson pg 59
k = 1./(pi*A*e);
P_wceil = @(W_S) (((Vv./(sqrt(((2/rhoceil).*(W_S)).*sqrt(k/(3*CDmin))))) + (4*sqrt((k*CDmin)/3)))*V_cr)/(550*np);
Psh_Wceil = P_wceil(W_S);


%% GRAPHING SHIT

%%% Graph %%%
figure(1)
% Takeoff Distance
plot(Range,WP(Range),'m'); hatchedline(Range,WP(Range),'m',(300-180)*pi/180); hold on;
% plot(Range,WP1(Range),'m'); hatchedline(Range,WP1(Range),'m',(300-180)*pi/180); hold on;
WeightRatio=.84; % for turboprops

% LANDING
% 
plot(WS*ones(1,100),linspace(0,1.4,100)); hatchedline(WS*ones(1,100),linspace(0,1.4,100),'r-'); hold on;
% Specify line style, angle, aspect ratio, spacing, and length
linespec = 'r-';       % Red solid line
theta = 45  * pi / 180; % 45-degree hatch angle
ar = 1;                % Aspect ratio
spc = -0.02;           % Hatch spacing (fraction of x-range)
len = 1.2;             % Shorter hatch length relative to spacing

% Speed
plot(W_S,WshOverW_P); hatchedline(W_S,WshOverW_P,'b',theta,ar,-0.015,.02); hold on;
% Ceiling
plot(W_S,Psh_Wceil,'k'); hatchedline(W_S,Psh_Wceil,'k',theta,ar,-0.015,.02);hold on

%% PLot competitors
% TurboProp
for i =1: height(planeData)
    % Quick Conversions
    W=MTOW(i);
    WS=W/S(i);
    PW=P(i)/W;
    plot(WS,PW,'*'); hold on;
    legend(names)
end


% Graph Parameters
xlim([10,100]); xlabel("W / S"); ylabel("P/ W"); grid on; ylim([0,1.4]);
text(77,.3,'Takeoff Distance'); text(10,.15,'Ceiling'); 
text(77,1.2,"Landing Distance"); text(13,.6,'Min Cruise Speed');
 %ylim([0,5.0]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% %% %% %% %% %% %% %% %% %% JETS %% %% %% %% %% %% %% %% %% %% %% %%
AR_jet=5.15; % average of comp asses
Cl_maxL=2.5;
Cl_maxTO=2;

%% Takeoff
% TOP_jet=WS/(sigma*Cl_maxTO*(TW));
S_TOFL=2350; % [ft] - from comp assesment
TOPNUM=37.5;
TOP=S_TOFL/ TOPNUM;

TW_TO= @(WS) WS/(sigma*Cl_maxTO*TOP);
figure(2)
plot(Range,TW_TO(Range),'m'); hatchedline(Range,TW_TO(Range),'m'); hold on;

%% Landing
S_L=2670; % [ft] - from comp assesment
S_FL=S_L/.6;
V_a=sqrt(S_FL/.3); 
V_s=V_a/1.3;
V_s=knots2ftperS(V_s);
WS=((V_s^2 *rho * Cl_maxL)/2)/.84;
plot(WS*ones(1,100),linspace(0,1.4,100),'r'); hatchedline(WS*ones(1,100),linspace(0,1.4,100),'r-'); hold on;

%% Cruise Speed
%cruise speed at sea level
T_Wcr= @(W_S) ((C_D_0*q./(W_S)+((W_S)*1./(q*pi*A*e))));
T_Wcruise = T_Wcr(W_S);
%ceiling at 25000 ft
Vv =  1.667;% ft/s obtiained from gudmunson no clue were they got pg 59
CDmin = .025; %value obtained fromchart in gudmunson pg 59
plot(WS, T_Wcruise,'b'); hatchedline(W_S, T_Wcruise, 'b',theta,ar,-0.015,.02); hold on
%% Ceiling
k = 1./(pi*A*e);
T_Wce = @(W_S) (Vv./(sqrt(((2/rhoceil).*(W_S)).*sqrt(k/(3*CDmin))))) + (4*sqrt((k*CDmin)/3));
T_Wceiling = T_Wce(W_S);
plot(W_S,T_Wceiling,'k'); hatchedline(W_S, T_Wceiling, 'k',theta,ar,-0.015,.02);

%% PLot competitors
% TurboJet
for i =1: height(planeData)
    % Quick Conversions
    W=MTOWJ(i);
    WS=W/SJ(i);
    TW=Thrust(i)/W;
    plot(WS,TW,'*'); hold on;
end

grid on;
xlim([10,100]); xlabel("W / S"); ylabel("T/ W"); grid on; ylim([0,1.4]);
text(45,.05,'Takeoff Distance'); text(80,.37,'Ceiling'); 
text(65,1.25,"Landing Distance"); text(13,.6,'Min Cruise Speed');




