% LANDING DISTANCE


% STEP 1 - TOTAL LANDING DISTANCE
S_L=1500;
landingRelationship= 1.938; %% ASSUMPTION FROM ROSKAM FROM DATA %%
S_LG=S_L/landingRelationship; % [knots]


% STEP 2 - LANDING ROLL RELATIONSHIP

RoskamSpeed2Dist=.265; % Pretty closely verified through our data

V_S=sqrt(S_LG/RoskamSpeed2Dist); % knots 
V_a= knots2ftperS(1.3 * V_S); % Calc V_S & convert knots to ft/s


Range=linspace(1,100,100); % range of W/P values


Cl_maxL= @(WS) (2*WS)/(rho*V_a^2)



%%% Landing %%%
WeightRatio=.84; % for turboprops
% WS= (rho*Cl_maxL(WS)*V_a^2)/2
Cl_maxL=2;
WS1= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio
Cl_maxL=2.2;
WS2= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio
Cl_maxL=2.6;
WS3= ((rho*Cl_maxL*V_a^2)/2)/WeightRatio

