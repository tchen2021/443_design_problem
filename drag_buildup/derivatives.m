%% Longitudinal Aerodybamic Deviatives Calculator 
%Tyson Chen


%takes geometric data of the wing and fuel as specified in drag buildup
%notes on Canvas and calculates 

%Wing: lift slope (alpha_w), zero lift angle of attack (alpha_0W),
%down-wash derviatives (delepsilon/delalpha)

%Fuselage/Nacelles: lift slope (alpha_B), zero lift angle of attack
%(alpha_0B)

%Tail: lift slope (alpha_T = alpha_1), aero lift angle of attack (alpha_0T)

%houskeeping
clc; clear; close all;

%% inputting wing geometric parameters
%wing
S = 121;               %area [ft^2]
cbar = ;            %mean geometric chord [ft]
cbarbar = ;         %mean aerodynamic chord [ft]
xbar = ;            %position of m.a.c [ft]
ybar = ;            %position of m.a.c [ft]
Lambda_0 = ;          %sweep angle at the l.e 
Lambda_quarter = ;   %sweep angle at 1/4c
Lambda_half = ;      %sweep angle at 1/2c
Lambda_te = ;        %sweep angle at t.e
%control surfaces
S_aileron = ;               %area of aileron [ft^2]
cbar_aileron = ;            %mean geometric chord [ft]
cbarbar_aileron = ;         %mean aerodynamic chord [ft]
S_elevator = ;              %area of elevator [ft^2]
cbar_elevator = ;           %mean geometric chord [ft]
cbarbar_elevator = ;        %mean aerodynamic chord [ft]
S_rudder = ;                %area of rudder [ft^2]
cbar_rudder = ;             %mean geometric chord [ft]
cbarbar_rudder = ;          %mean aerodynamic chord [ft]
i_T = ;                     %tail incidence 

%fuselage
lf = 20.333;              %length of fuselage [ft]
la = ;              %length of nose [ft]
lb = ;              %length of cylindrical portion [ft]
lc = ;              %length of tail [ft]

D = ;               %max diameter of the EBR


%need to calculate max diameter, location of max diameter, wetted area of
%nose, cylindrical portion, and the tail

%% lift slope
M = ;                       %Mach number
phiprime_TE = ;             %typical airfoil trailing edge angle
K = ;                       %from interpolation of fig. B.1.1(a)
C_l_alpha_theory = ;        %interpolate from fig. B.1.1(b)
beta = sqrt(1-M.^2);        %sideslip angle 
C_l_alpha = (1.05./beta).*K.*C_l_alpha_theory;  %
kappa = beta*C_l_alpha/(2*pi);
C_L_alpha = ;               %interpolate from fig. B.1.2

%% zero-lift angle of attack
alpha_0W_root = ; %zero lift angle of attack, 2D, at the wing root
iprime_r = ;      %geometric incidence of the wing at the root
J = ;             %empirical factor
epsilon = ;       %aerodynamic twist of the wing
lambda = ;        %taper ratio c_t/c_s
alpha_0W = alpha_0W_root - iprime_r + J*epsilon; %zero-lift angle of attack

%% Downwash
K_A = ;     %wing-aspect-ratio factor obtained from fig B.5.1 
K_lambda = ;%wing-taper-ratio factor obtained from fig B.5.2
K_H = ;     %horizontal-tail-location factor from fig B.5.3
delepsilondelalpha = 4.44.*(K_A.*K_lambda.*K_H.*cos(Lambda_quarter).^0.5).^1.19;

%% Fuselage Lift
alpha = ;       %angle of attack of the airplane 
alpha_0B = ;    %zero lift angle of fuselage
K_fuse = ;      %apparent mass function, interpolated from figure
eta =    ;      %relation between drag of an infinity cylinder and finite cylinder
C_DC = ;        % drag coefficient for inclined flow over an infinity cylinder
S_p_x0 = ;      %planform area of the EBR behind a point x0, definied as point where fuselage starts to have shedding vortex
C_LB = *(alpha-alpha_0B) ./ S) .* ( ((K.*pi.*D.^2) / (2)) + eta.*C_DC .*(alpha-alpha_0B).*S_p_x0 ); %lift of a fuse

%wing body lift
k_wb = ;        %factor of interference of wing-body for lift coefficient obtained from fig
k_bw = ;        %factor of interference of wing-body for lift coefficient obtained from fig
Se = ;          %wing area not covered by fuselage [ft^2]
C_LWB = C_LB + (k_wb-k_bw).*a_w.*(alpha-alpha_0W).*(Se./S);     %total lfit of wing-body system, alpha in rad

%horizontal tail lift
a1 = ;          %lift slope of the tail
C_LT = a1 .* (alpha.*(1-delepsilondelalpha) - (epsilon_0 + i_t)) .* (S_T/S);

%total lift
C_L = C_LW + C_LT;

%% Wing Drag
V = 

Re = (V .* rho .* MAC) ./ mu;                   %calculate Reynold's number depending on flight conditions (speed,alittude, etc)

alpha_local = alpha + iprime_r + epsilon_i;     %calculate local angle of attack of each section

%given CL calculated above and Re, CD can be obtained using airfoil data
%from Theory of Wing Sections or airfoiltoolbox

%%%%%%%%%%%%%%%
%organize and export a table of CD depending on section, speed, weight, CG
%position, configuration, altitude, assuming elliptical lift distribution

%%%%%%%%%%%%%%%
C_D0W = 2.* ((sum(C_Di.*S_i)) /(sum(S_i))) * (S_e/S);   %total drag coefficient

delta1 = ;          %interpolated from figure
delta2 = ;          %interpolated from figure
C_DiW = ((C_LW.^2) ./ (pi.*S)) .* (1+delta1+delta2);    %induced drag coefficient of wing


%Tail
iprime_T = ;        %geometric incidence of the horizontal tail
alpha_T = alpha + iprime_T - delepsilondelalpha .* alpha;   %the angle of attack of horizontal tail
alpha_V = 0;        %alpha on the vertical tail is 0 for straight leveled flight
C_D0T = 2.* ((sum(C_DTi.*S_Ti)) /(sum(S_Ti))) * (S_Te/S);   %subscript T refers to the tail (horizontal or vertical)

delta3_f1 = ;          %interpolated from figure flap 1
delta3_f2 = ;          %interpolated from figure flat 2

delta_C_D0T = delta1 * delta2 * (delta3_f1 - delta3_2);
