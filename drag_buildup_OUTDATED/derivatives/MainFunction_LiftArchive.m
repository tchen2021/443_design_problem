%% Longitudinal Aerodybamic Deviatives Calculator 
%Tyson Chen, Roman Niemiec


%takes geometric data of the wing and fuel as specified in drag buildup
%notes on Canvas and calculates 

%Wing: lift slope (alpha_w), zero lift angle of attack (alpha_0W),
%down-wash derviatives (delepsilon/delalpha)

%Fuselage/Nacelles: lift slope (alpha_B), zero lift angle of attack
%(alpha_0B)

%Tail: lift slope (alpha_T = alpha_1), aero lift angle of attack (alpha_0T)

%houskeeping
clc; clear; close all;

%adding derivatives functions and data
addpath derivatives;
%% importing flight conditions
h = 3e4;                    %to be parameterized
V = 280 * 1.68781;          %knots to ft/s
[~, T, rho, a] = atmosphere(h);                    %30,000 ft cruise
M = V/a;                       %Mach number

%% inputting wing geometric parameters
%wing
S = 190;                %area [ft^2]
cbar = 5.45;            %mean geometric chord [ft]
C_tip = 3.5;            %wingtip chord [ft]
C_root = 7.4;           %root chord [ft]
lambda = C_tip/C_root;  %taper ratio
cbarbar = (((2/3) * ((C_root+C_tip-((C_root*C_tip)/(C_root+C_tip))) / (1 + lambda))) + cbar)/2; %mean aerodynamic chord [ft]
%xbar = ;                 %position of m.a.c [ft]
%ybar = ;                 %position of m.a.c [ft]
Lambda_0 = 8.776;         %sweep angle at the l.e [deg]
Lambda_quarter = 8.776;   %sweep angle at 1/4c
Lambda_half = 8.776;      %sweep angle at 1/2c
Lambda_te = 8.776;        %sweep angle at t.e
AR = 6.25;
b = 35.35;                %wing span [ft]
l_H = 14.59404;           %wing m.a.c to horizontal tail m.a.c [ft]
t_root = 0.15 * cbar;                %thickness of the root, assumed 0.15c 
% t_tail = 0.10 * ;                 %thickness of the tail           
S_e = S * (1-AR^-1);            %effective wing area that contributes to downwash
%control surfaces
S_aileron = 0.25 * 0.225 * S;               %area of aileron [ft^2] rough estimation b_a ~ 0.25b c_a ~ 0.225c
cbar_aileron = 0.225 * cbar;                %mean geometric chord [ft]
cbarbar_aileron = 0.225 * cbarbar;          %assume to taper mean aerodynamic chord [ft]
S_elevator = 22.93844;      %area of elevator [ft^2]
cbar_elevator = (38.39/12) * 0.25;           %mean geometric chord [ft]
cbarbar_elevator = cbar_elevator;            % assume no taper mean aerodynamic chord [ft]
S_rudder = 19.69203;        %area of rudder [ft^2]
cbar_rudder = (48.58/12) * 0.25;             %mean geometric chord [ft]
cbarbar_rudder = cbar_rudder;          %mean aerodynamic chord [ft]
i_T = 0;  %tail incidence 

%wing-tail
h_H = 31.5/12;      %vertical distance between wing mac and tail mac
%fuselage
lf = 20.333;              %length of fuselage [ft]
la = 7.3025;              %length of nose [ft]
lb = 13.2958;             %length of cylindrical portion [ft]
lc = 14.1866;             %length of tail [ft]

%% fuselage EBR

%D =;               %max diameter of the EBR
S_A_wet = 27.09;    %wetted area of the nose area [ft^2]
S_B_wet = 101.23;   %wetted area of the center area [ft^2]
S_C_wet = 53.73;    %wetted area of the tail area [ft]^2]

%need to calculate max diameter, location of max diameter (optional)

% Import data file
DragBuildUp = importfileTEST('DragBuildUp.xlsx');
% Take individual variables from excel sheet
for i=1:height(DragBuildUp)
    assignin('base',DragBuildUp{i,2},DragBuildUp{i,3});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          LIFT SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift slope - Wing
% Solve Graph B11A for K value
[K_Wing] = B11A(h,V,C_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Wing] = B11B(t_root,t_tail,C_root,C_tip);

%%% 5. Calculate beta
[M_Wing] = machCalc(V,h);
beta_Wing=sqrt(1-M_Wing^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Wing=(1.05/beta_Wing)*K_Wing*Cl_alpha_theory_Wing;

%%% 6. Calculate k
k_Wing=(beta_Wing*Cl_alpha_Wing)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);
% [CL_alpha_Wing] = B12(xaxis)/A;
Cl_alpha_Wing= ( (2*pi)/(2+sqrt( ((A^2 * beta_Wing^2)/2)*(1+ (tand(Lambda_half)^2)/beta_Wing^2 )+4))) /A; % Function mentioned on graph 

%% Lift slope - TAIL SECTION %%
% Solve Graph B11A for K value
[K_Horiz] = B11A(h,V,C_Vert_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Horiz] = B11B(t_Vert_root,t_Vert_tip,C_Vert_root,C_Vert_tip);

%%% 5. Calculate beta
[M_Horiz] = machCalc(V,h);
beta_Horiz=sqrt(1-M_Horiz^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Horiz=(1.05/beta_Horiz)*K_Horiz*Cl_alpha_theory_Horiz;

%%% 6. Calculate k
k_Horiz=(beta_Horiz*Cl_alpha_Horiz)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);
% [CL_alpha_Wing] = B12(xaxis)/A;
Cl_alpha_Horiz= ( (2*pi)/(2+sqrt( ((A^2 * beta_Horiz^2)/2)*(1+ (tand(Lambda_Tail_half)^2)/beta_Horiz^2 )+4))) /A; % Function mentioned on graph 



%% zero-lift angle of attack
alpha_0W_root = deg2rad(-3); %zero lift angle of attack, 2D, at the wing root
iprime_r = deg2rad(2);      %geometric incidence of the wing at the root
%J = ;             %empirical factor
epsilon = 0;       %aerodynamic twist of the wing
alpha_0W = zeroliftalpha(iprime_r, epsilon, alpha_0W_root, AR, lambda);
%alpha_0W = alpha_0W_root - iprime_r + J*epsilon; %zero-lift angle of attack

%% Downwash
armfactor = 2*l_H / b; % Ratio of mac distance over area
tailfactor = abs(2*h_H/b);
delepsilondelalpha = downwash(Lambda_quarter, AR, lambda, tailfactor, armfactor,h_H,b,l_H); % how the wake effects of the wing affect the flow of the tail


%% Fuselage Lift
%largest fuse area 9786.86
% 2786.86
A = 2786.86 / 144;
D = sqrt(4*A/pi);
S_P_x0 = 934.73 / 144;
alpha = 5;       %angle of attack of the airplane 
alpha_0B = 0;    %zero lift angle of fuselage
lf_D = lf / D;


C_LB = fuselift(lf_D, alpha, alpha_0B, S, D, S_P_x0, M, la, lf);

%%% Wing body lift %%%
Se_S = 0.8;     %ratio of wing not covered by the fuselage to the total wing area, estimation from online
D_b = D/b;
alpha_W = 1;
C_LW = wingbodylift(C_LB, alpha_W, alpha, alpha_0W, Se_S, D_b);
% k_wb = ;        %factor of interference of wing-body for lift coefficient obtained from fig
% k_bw = ;        %factor of interference of wing-body for lift coefficient obtained from fig
% Se = ;          %wing area not covered by fuselage [ft^2]
% C_LWB = C_LB + (k_wb-k_bw).*a_w.*(alpha-alpha_0W).*(Se./S);     %total lfit of wing-body system, alpha in rad

%%% horizontal tail lift %%%
a1 = Cl_alpha_Horiz;          %lift slope of the tail
epsilon_0=0;
C_LT = a1 .* ((alpha*pi/180).*(1-delepsilondelalpha) - (epsilon_0 + i_T)) .* (S_Vert/S); % Calculate with alpha in rads

%total lift
C_L = C_LW + C_LT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          DRAG SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wing Drag
% V = 
% 
% Re = (V .* rho .* MAC) ./ mu;                   %calculate Reynold's number depending on flight conditions (speed,alittude, etc)
% 
% alpha_local = alpha + iprime_r + epsilon_i;     %calculate local angle of attack of each section
% 
% %given CL calculated above and Re, CD can be obtained using airfoil data
% %from Theory of Wing Sections or airfoiltoolbox
% 
% %%%%%%%%%%%%%%%
% %organize and export a table of CD depending on section, speed, weight, CG
% %position, configuration, altitude, assuming elliptical lift distribution
% 
% %%%%%%%%%%%%%%%
% C_D0W = 2.* ((sum(C_Di.*S_i)) /(sum(S_i))) * (S_e/S);   %total drag coefficient

C_D0W = iterator(airfoil,alpha,iprime_r, epsilon,cbarbar, S_arr,V,rho,mu,Se, S);
C_DiW = C_DiWfun(lambda, Lambda_quarter, C_LWB-C_LB, AR);


%Tail
iprime_T = 0;        %geometric incidence of the horizontal tail assume 0
alpha_T = alpha + iprime_T - delepsilondelalpha .* alpha;   %the angle of attack of horizontal tail
alpha_V = 0;        %alpha on the vertical tail is 0 for straight leveled flight
C_D0T = 2.* ((sum(C_DTi .* S_Ti)) / (sum(S_Ti))) * (S_Te/S);   %subscript T refers to the tail (horizontal or vertical)

% delta3_f1 = ;          %interpolated from figure flap 1
% delta3_f2 = ;          %interpolated from figure flat 2
% 
% delta_C_D0T = delta1 * delta2 * (delta3_f1 - delta3_2); %additional drag caused by a plain flap (elevator)


C_DiT = C_DiTfun(C_LT, A_T, S, S_T, t_c, beta, c_C);


%% Fuselage and Nacelles
C_D0F = C_D0Ffun(Re, S, S_A_wet, S_B_wet, S_C_wet, S_CAB, 't', l_f, D, Lambda_f);       %assumes smooth surface turbulent flow
C_DiB = C_LB*(alpha-alpha_0B);          %induced drag coefficient of the fuselage w/ alpha in radians

%% Interference Drag

%wing-fuselage
deltaC_DWf = 0.05 * (C_D0F + C_DiF); 

%wing-nacelle
n = ;   %number of nacelles
deltaC_Dwn = 0.05 * (C_D0N + C_DiN) * n;

%Tail-Fuselage
n_c = ;     %number of corners between tail and fuselage
c_j = ;     %chord of the tail at the junction
deltaC_Dtf = n_c * (0.8 * (t/c)^3 - 0.0005) * (c_j^2/S);

deltaC_Dtt = (n_c/2) * (17*(t/c)^4 - 0.05*(t/c)^2) * (c_j^2/S);

%cooling and air intake
mdot = ;        %mass flow rate of air inside the intake
q = ;           %dynamic pressure
deltaV =  ;     %variation of the velocity between the intake and the outlet
C_D_cooling = (mdot*deltaV)/(q*S);


