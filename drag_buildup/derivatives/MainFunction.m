%% Longitudinal Aerodybamic Deviatives Calculator 
%Authors: Tyson Chen & Roman Niemiec
%This main script takes geometric data of the plane design pulled from an excel sheet named
% "DragBuildUp.xls and calculates aerodynamic peformance of it using the drag buildup
% notes on Canvas
clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                          DATA SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding derivatives functions and data
addpath derivatives;
% Import data file
DragBuildUp = importfileTEST('DragBuildUp.xlsx');
% Take individual variables from excel sheet
for i=1:height(DragBuildUp)
    assignin('base',DragBuildUp{i,2},DragBuildUp{i,3});
end
                    
j=1;
airfoil='NACA 2412';
%%% ITERATOR for CL CD 
for alpha=begin:ending 

% Get variables from Speed and altitude
[p, T, rho, a] = atmosphere(h); % get temp in [R]
[M] = machCalc(V,h); % Mach number based on speed of air outside airplane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          LIFT SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift slope - Wing
% Solve Graph B11A for K value
[K_Wing,Re] = B11A(h,V,C_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Wing] = B11B(t_root,t_tail,C_root,C_tip);

%%% 5. Calculate beta
beta_Wing=sqrt(1-M^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Wing=(1.05/beta_Wing)*K_Wing*Cl_alpha_theory_Wing;

%%% 6. Calculate k
k_Wing=(beta_Wing*Cl_alpha_Wing)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);
% [CL_alpha_Wing] = B12(xaxis)/A;
% CL_alpha_Wing= ( (2*pi)/(2+sqrt( ((A^2 * beta_Wing^2)/2)*(1+ (tand(Lambda_half)^2)/beta_Wing^2 )+4))) *A; % Function mentioned on graph 
CL_alpha_Wing= (A*2*pi) / (2+ sqrt( ((A^2)*(beta_Wing^2))/k_Wing^2 *(1+(tand(Lambda_half)^2)/beta_Wing^2)+4))

%% Lift slope - TAIL SECTION %%
% Solve Graph B11A for K value
[K_Horiz,Re] = B11A(h,V,C_Vert_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Horiz] = B11B(t_Vert_root,t_Vert_tip,C_Vert_root,C_Vert_tip);

%%% 5. Calculate beta
beta_Horiz=sqrt(1-M^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Horiz=(1.05/beta_Horiz)*K_Horiz*Cl_alpha_theory_Horiz;

%%% 6. Calculate k
k_Horiz=(beta_Horiz*Cl_alpha_Horiz)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);
% [CL_alpha_Wing] = B12(xaxis)/A;
CL_alpha_Horiz= ( (2*pi)/(2+sqrt( ((A^2 * beta_Horiz^2)/2)*(1+ (tand(Lambda_Tail_half)^2)/beta_Horiz^2 )+4))) *A; % Function mentioned on graph 

%% zero-lift angle of attack
% alpha_0W_root = deg2rad(-3); %zero lift angle of attack, 2D, at the wing root
% iprime_r = deg2rad(-5);      %geometric incidence of the wing at the root
%J = ;             %empirical factor
epsilon = 0;       %aerodynamic twist of the wing
alpha_0W = zeroliftalpha(iprime_r, epsilon, alpha_0W_root, AR, lambda);

%% Downwash
armfactor = 2*l_H / b; % Ratio of mac distance over area
tailfactor = abs(2*h_H/b);
delepsilondelalpha = downwash(Lambda_quarter, AR, lambda, tailfactor, armfactor,h_H,b,l_H); % how the wake effects of the wing affect the flow of the tail

%% Fuselage Lift
%largest fuse area 9786.86
% 2786.86
A_area = 2786.86 / 144;
D = sqrt(4*A_area/pi);
S_P_x0 = 934.73 / 144;

alpha_0B = 0;    %zero lift angle of fuselage
lf_D = lf / D;

% Solve Lift value for lift body
C_LB = fuselift(lf_D, alpha, alpha_0B, S, D, S_P_x0, M, la, lf);

%%% Wing body lift %%%
Se_S = 0.8;     %ratio of wing not covered by the fuselage to the total wing area, estimation from online
D_b = D/b;
a_W = CL_alpha_Wing;
C_LW = wingbodylift(C_LB, a_W, alpha, alpha_0W, Se_S, D_b); %alpha in rad

%%% horizontal tail lift %%%
a1 = CL_alpha_Horiz;          %lift slope of the tail
epsilon_0=0;
C_LT = a1 .* ((deg2rad(alpha)*(1-delepsilondelalpha)) - (epsilon_0 + i_T))* (S_T/S); % Calculate with alpha in rads

%total lift
C_L(j) = C_LW + C_LT;
% end % end of test for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          DRAG SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wing Drag
% For loop to iterate through 4 sections of our wing
S_i=[S1,S2,S3,S4];
cbarbari=[cbarbar1,cbarbar2,cbarbar3,cbarbar4];
% 
% for i =1; i<=4; 
% %     [mu] = SutherlandsEquation(h);
% %     %calculate local Reynold's number
% %     Re_local(i) = (V* rho * cbarbari(i) )/ mu;                 
% %     %calculate local angle of attack of each section
% %     alpha_local(i) = alpha + iprime_r + epsilon_i;
% 
% %     % Local Drag
% %     C_Di= ;
%      
% end % end of for loop

C_D0W = iterator(airfoil,alpha,iprime_r, epsilon,cbarbari, S_i,V,rho,Se, S,h);
C_DiW = C_DiWfun(lambda, Lambda_quarter, C_LB-C_LB, AR);




%Tail
iprime_T = 0;        %geometric incidence of the horizontal tail assume 0
alpha_T = alpha + iprime_T - delepsilondelalpha .* alpha;   %the angle of attack of horizontal tail
alpha_V = 0;        %alpha on the vertical tail is 0 for straight leveled flight
% C_D0T = 2.* ((sum(C_DTi .* S_Ti)) / (sum(S_Ti))) * (S_Te/S);   %subscript T refers to the tail (horizontal or vertical)
C_D0T=0; %  WE AINT DOING THIS YET
% delta3_f1 = ;          %interpolated from figure flap 1
% delta3_f2 = ;          %interpolated from figure flat 2
% 
% delta_C_D0T = delta1 * delta2 * (delta3_f1 - delta3_2); %additional drag caused by a plain flap (elevator)


C_DiT = C_DiTfun(C_LT, A_T, S, S_T, t_c, beta_Horiz, cf_C); % UNCOMMENT SOON


%% Fuselage and Nacelles
C_D0F = C_D0Ffun(Re, S, S_A_wet, S_B_wet, S_C_wet, S_CAB, 't', lf, D, Lambda_f);       %assumes smooth surface turbulent flow
C_DiB = C_LB*(deg2rad(alpha)-(deg2rad(alpha_0B)));          %induced drag coefficient of the fuselage w/ alpha in radians

C_D(j)=C_DiB+C_DiT+C_D0T+C_DiW+C_D0W + C_D0F;

% % %% Interference Drag
% % 
% % %wing-fuselage
% % deltaC_DWf = 0.05 * (C_D0F + C_DiF); 
% % 
% % %wing-nacelle
% % n = ;   %number of nacelles
% % deltaC_Dwn = 0.05 * (C_D0N + C_DiN) * n;
% % 
% % %Tail-Fuselage
% % n_c = ;     %number of corners between tail and fuselage
% % c_j = ;     %chord of the tail at the junction
% % deltaC_Dtf = n_c * (0.8 * (t/c)^3 - 0.0005) * (c_j^2/S);
% % 
% % deltaC_Dtt = (n_c/2) * (17*(t/c)^4 - 0.05*(t/c)^2) * (c_j^2/S);
% % 
% % %cooling and air intake
% % mdot = ;        %mass flow rate of air inside the intake
% % q = ;           %dynamic pressure
% % deltaV =  ;     %variation of the velocity between the intake and the outlet
% % C_D_cooling = (mdot*deltaV)/(q*S);
% % 
j=j+1;
end % end of for loop through


% Check the data
% Plot CL vs Alpha
plot(begin:ending,C_L); xlabel("Alpha");ylabel("CL");grid on;