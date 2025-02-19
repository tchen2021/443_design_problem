% Roman Niemiec - Drag Build Up 
% WING SECTION %
% 1/15/2025

%% 1. Calculate typical airfoil trailing edge angle 
% Geometry
c=6; % Chord [ft]
Y90 = 1; % height of LE [ft]
R=10000000:10000000


% EQUATION B11-1 [TRAILING EDGE ANGLE]
theta_TE= atand(2*((Y90/2) / (.09*c)));

.5*tand(theta_TE);
Reynolds=log(R);


%% Zero-Lift Angle of Attack
% Initialize Variables
alpha_0W=1; %zero lift angle of attack of the wing
alpha_0W_Root=1; % zero-lift angle of attack, two dimensional, at the wing root. 
i_r=1; % geometric incidence of the wing at the root
J =1; % Empirical factor
Epsilon= 1; %aerodynamic twist of the wing
lambda=1; %c_t/c_s= lambda - taper ratio


%The zero lift angle of attack of a wing can be determined by:
alpha_0W= alpha_0W_Root - i_r + J*Epsilon;

% For a wing with flaps, the zero lift angle of attack will be change by
delta_alpha_0W = -delta_CL/ CL_alpha; 

%% Downwash %%
% Initialize Variables




% The average lowspeed downwash gradient at the horizontal tail is given by:
downWashGradient= 4.44*(K_A*K_epsilon*K_H*(cos()))
%% Part 4 - Lift %%
% [] = DragBuildUp_4_Lift()