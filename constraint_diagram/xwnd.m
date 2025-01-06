% Constraint Diagram - Group 5 - Roman Niemiec
clc;clear;close all;format short;
% Added friction values to takeoff calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint Diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pick which loading Method you desire
    % load planeData
    fileLocation="C:\Users\olive\Documents\Aero\Aero433\planeData.xlsx";
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
W_S =linspace(1, 100, 100); % lb/ft^2 wing loading
C_D_0 = 0.025;%assuming this for all of them can change

%%% Standard Atmo Values %%%
[p0, T0, rho0, a0] = atmosphere(0); % Calculate standard atmosphere values at sea level
[p, T, rho, a] = atmosphere(h);% Calculate standard atmosphere values at variable height
[~, ~, rhospeed, ~] = atmosphere(14000);% Calculate standard atmosphere values at variable height

% Correct for non standard atmosphere
sigma=rho./rho0;

%Graph setup
% Specify line style, angle, aspect ratio, spacing, and length
linespec = 'r-';       % Red solid line
theta = 45  * pi / 180; % 45-degree hatch angle
theta2 = -45  * pi / 180; % 45-degree hatch angle

ar = 1;                % Aspect ratio
spc = -0.02;           % Hatch spacing (fraction of x-range)
len = 1.2;             % Shorter hatch length relative to spacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Takeoff %%%
S_TOFL=2250; % no ground constant
Range=linspace(1,100,100); % range of W/P values
V=343; %calculating using clmax to and 22000%400; %ft/s
eff=.82;
Cl_maxTO=1.8; %%% ASSUMPTION
TOPchange= @(TOPtest) 8.134*TOPtest + .0149*TOPtest^2 - S_TOFL;
TOP=fzero(TOPchange,.1);

% Equation for W/P (Power Loading)
% WP= @(WS) WS/(sigma*Cl_maxTO*TOP); SUBJECT TO DESIGN REVIEW TEST

%%% FRICTION CALC VALUES %%%
ip = 5.75;%from roskam typical value for constant speed props which are more efficient
s_TOGround = 2250;         % Takeoff ground run distance
k1 =.0376;              % Constant k1
k2 =ip*((1/13)^(1/3)) ;    %proppeller disk loading range from 10-30  
mu_G = .3;     % Ground friction coefficient soft ground worst case
mu_G_med=.08;
mu_Gmin=.02;
% Calculate (X/W)_TO
% Best Friction case (concrete/asphalt)
WP = @(WS) ( (k1 * (WS)) / (s_TOGround * rho * Cl_maxTO) + mu_Gmin + 0.72 * C_D_0 ) / k2;
% % Medium best case (wet grass .08)
P_W_TOground_med = ( (k1 * (W_S)) / (s_TOGround * rho * Cl_maxTO) + mu_G_med + 0.72 * C_D_0 ) / k2;
% Worst Case (soft ground) 
P_W_TOground = ( (k1 * (W_S)) / (s_TOGround * rho * Cl_maxTO) + mu_G + 0.72 * C_D_0 ) / k2;
% using gudmunson formula, the answer was all "fucked up"
% T_Wgroudngud = (((V_LOF^2 ./ (2 * g * S_G)) + (qgud * C_D_TO ./ W_S) + mu * (1 - (qgud * CL_max_TO ./ W_S)))*V_LOF)/(550*np);

%%% Graph %%%
figure(1)
% Takeoff Distance
plot(Range,WP(Range),'b'); hatchedline(Range,WP(Range),'b',(300-180)*pi/180,ar,-0.015,.01); hold on;

plot(W_S,P_W_TOground_med,'b'); hatchedline(W_S,P_W_TOground_med,'b',(300-180)*pi/180,ar,-0.015,.01); hold on; % Friction for ground
%plot(W_S,P_W_TOground,'b');
%hatchedline(W_S,P_W_TOground,'b',(300-180)*pi/180); hold on; % Friction
%for ground wrst case

%% landing Calcs
WeightRatio=.84;% for turboprops
%1500 ft case
% Landing Distance
S_L=2350;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]
% STEP 2 - LANDING ROLL RELATIONSHIP
RoskamSpeed2Dist=.265; % Pretty closely verified through our data
V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s
Cl_maxL=2.4;
WS= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio;

%%% Graph %%% 
plot(WS*ones(1,100),linspace(0,1.4,100)); hatchedline(WS*ones(1,100),linspace(0,1.4,100),'r-',(300-180)*pi/180); hold on;


%competetitor
S_L=3000;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]
% STEP 2 - LANDING ROLL RELATIONSHIP
RoskamSpeed2Dist=.265; % Pretty closely verified through our data
V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s
Cl_maxL=2.4;
WS= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio;
%%% Graph %%% 
plot(WS*ones(1,100),linspace(0,1.4,100)); hatchedline(WS*ones(1,100),linspace(0,1.4,100),'r-',(300-180)*pi/180); hold on;




%% speed requirement 
% min cruise speed should be greater than 280 knotts
% values most likely needed q dynamic density at sea level, W_S wing loading, 
% assuming a AR, and e rn

%jets have not done yet
W_S =linspace(1, 100, 100); % lb/ft^2 wing loading

 V_cr = 472.587;%ft/s
 A = 6.6;%assumption from competitor - average of comp assesment
 e =0.8; % assumption for all of them
 q= (1/2)*rhospeed*V_cr^2;%same q for all of them
C_D_0 = 0.025;%assuming this for all of them can change
TOverW = @(W_S) (C_D_0*q./(W_S)+W_S*1./(q*pi*A*e));%jets
TOverW_Speed = TOverW(W_S);

% proppeler
np = 0.80;%michaels value

WshaftverP = @(W_S) (((C_D_0*q./(W_S)+(W_S)*1./(q*pi*A*e)))*V_cr)./(550*np);
WshOverW_P = WshaftverP(W_S);
%%% Graph %%%
plot(W_S,WshOverW_P,'g'); hatchedline(W_S,WshOverW_P,'g',theta,ar,-0.015,.02); hold on;



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
%%% Graph %%%
plot(W_S,Psh_Wceil,'k'); hatchedline(W_S,Psh_Wceil,'k',theta,ar,-0.015,.02);hold on

%% Plot wingloading at different P/W curve  porps

w_analysisprop=[2];

P_Wcurve = 1600/w_analysisprop;





%% PLot competitors
% TurboProp
for i =1: 2
    % Quick Conversions
    W=MTOW(i);
    WS=W/S(i);
    PW=P(i)/W;
    plot(WS,PW,'k*', 'MarkerSize', 9); hold on;
end

 

% Graph Parameters
xlim([30,100]); xlabel("W / S [ lb / ft^{2} ] ","FontSize",26); ylabel("P/ W [ Shp / lb ] ","FontSize",26); grid on; ylim([0,.3]);
% text(78,.3,'Takeoff Distance',"FontSize",28); text(82,.05,'Ceiling',"FontSize",28); 
% text(78,.78,"Landing Distance","FontSize",28); text(21,.6,'Min Cruise Speed',"FontSize",28);
set(gca, 'FontSize', 26);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% %% %% %% %% %% %% %% %% %% JETS %% %% %% %% %% %% %% %% %% %% %% %%
AR_jet=7.0; % average of comp asses
Cl_maxL=2.5;%2.4
Cl_maxTO=2.0;%1.8

%% Takeoff
% TOP_jet=WS/(sigma*Cl_maxTO*(TW));
S_TOFL=2250; % [ft] - from comp assesment no ground constant 
TOPNUM=37.5;
TOP=S_TOFL/ TOPNUM;

TW_TO= @(WS) WS/(sigma*Cl_maxTO*TOP);


% takoff using ground friction for jets
% ip = 5.75 for sonstant speed prop, 4.60 for fixed pitch props
% Given values (replace these with actual values)
s_TOGround = 2250;         % Takeoff ground run distance with grounf constant
k1 =.0447;              % Constant k1
lambda = 1.5;%engine bypass ratio
k2 = 0.75*((5+lambda)/(4+lambda)) ;    %proppeller disk loading range from 10-30  
mu_G = .08;     % Ground friction coefficient soft ground worst case
mu_G2 = .02;     % Ground friction coefficient soft ground worst case

% Calculate (X/W)_TO
T_W_TOground8 = ( (k1 * (W_S)) / (s_TOGround * rho * Cl_maxTO) + mu_G + 0.72 * C_D_0 ) / k2;

T_W_TOground2 = ( (k1 * (W_S)) / (s_TOGround * rho * Cl_maxTO) + mu_G2 + 0.72 * C_D_0 ) / k2;

figure(2)

%jet non ground
plot(Range,T_W_TOground2,'b'); hatchedline(Range,T_W_TOground2,'b',theta,ar,-0.015,.02); hold on;
plot(Range,T_W_TOground8,'b' ); hatchedline(Range,T_W_TOground8,'b',theta,ar,-0.015,.02);hold on;
%ground

% takoff using ground friction for props
% % ip = 5.75 for sonstant speed prop, 4.60 for fixed pitch props
% % Given values (replace these with actual values)
% ip = 4.60;%power ,oading
% s_TOG = 2250;         % Takeoff ground run distance
% k1 =.0376;              % Constant k1
% k2 =ip*((1/20)^(1/3)) ;    %proppeller disk loading range from 10-30  
% Cl_maxTO = 2.0; % Maximum lift coefficient at takeoff
% mu_G = 0;     % Ground friction coefficient soft ground worst case
% 
% % Calculate (X/W)_TO
% P_W_TOground = ( (k1 * (W_S)) / (s_TOG * rho * Cl_maxTO) + mu_G + 0.72 * C_D_0 ) / k2;

%% Landing
S_L=2350;%2700; % [ft] - from comp assesment
S_FL=S_L/.6;
V_a=sqrt(S_FL/.3); 
V_s=V_a/1.3;
V_s=knots2ftperS(V_s);
WS=((V_s^2 *rho * Cl_maxL)/2)/.84;
plot(WS*ones(1,100),linspace(0,1.4,100),'r'); hatchedline(WS*ones(1,100),linspace(0,1.4,100),'r-',(300-180)*pi/180); hold on;

%competeitilanddist
S_L=2900;%2700; % [ft] - from comp assesment
S_FL=S_L/.6;
V_a=sqrt(S_FL/.3); 
V_s=V_a/1.3;
V_s=knots2ftperS(V_s);
WS=((V_s^2 *rho * Cl_maxL)/2)/.84;
plot(WS*ones(1,100),linspace(0,1.4,100),'r'); hatchedline(WS*ones(1,100),linspace(0,1.4,100),'r-',(300-180)*pi/180); hold on;


%% Cruise Speed
%cruise speed at sea level
T_Wcr= @(W_S) ((C_D_0*q./(W_S)+((W_S)*1./(q*pi*AR_jet*e))));
T_Wcruise = T_Wcr(W_S);
%ceiling at 25000 ft
Vv =  1.667;% ft/s obtiained from gudmunson no clue were they got pg 59
CDmin = .025; %value obtained fromchart in gudmunson pg 59
plot(WS, T_Wcruise,'g'); hatchedline(W_S, T_Wcruise, 'g',theta,ar,-0.015,.02); hold on
%% Ceiling
k = 1./(pi*A*e);
T_Wce = @(W_S) (Vv./(sqrt(((2/rhoceil).*(W_S)).*sqrt(k/(3*CDmin))))) + (4*sqrt((k*CDmin)/3));
T_Wceiling = T_Wce(W_S);

plot(W_S,T_Wceiling,'k'); hatchedline(W_S, T_Wceiling, 'k',theta,ar,-0.015,.02);

%% Plot wingloading at different T/W curve jet

w_analysisjet=[2];

T_Wcurve = 1600/w_analysisjet;


%% PLot competitors
% TurboJet
for i =1: height(planeData)
    % Quick Conversions
    W=MTOWJ(i);
    WS=W/SJ(i);
    TW=Thrust(i)/W;
    plot(WS,TW,'k*', 'MarkerSize', 9); hold on;
end

grid on;
xlim([15,100]); xlabel("W / S [ lb / ft^{2} ] ","FontSize",26); ylabel("T / W [ lbf / lb] ","FontSize",26); grid on; ylim([0,.6]);
% text(45,.05,'Takeoff Distance'); text(80,.37,'Ceiling'); 
% text(65,1.25,"Landing Distance"); text(13,.6,'Min Cruise Speed');
set(gca, 'FontSize', 26);



%% sizing
% 
% CHt = (SHt*LHt)/(SW*Cmac);
% 
% CVt = S











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

