function [] = DragBuildUp_4_Lift()
% This is part 4 for the drag build up series to calculate the drag build
% up.
%   This section covers the lift of the fuselage and wings

%% 4 - Lift
%%% Initialize Variables %%%
% alpha = angle of attack in relaition to alpha_0B
% alpha_0B = fuselage zero lift angle of attack [GUESSED VALUE UNTIL CFD]
% S = Wing area
% K = apparent mass function
% D = maximum diameter of the EBR
% eta = relation between the drag of an infinity cylinder and finite cyldiner
% C_DC = drag coefficient for inclided flow over an infinity cylinder
% S_P_X0 = plan form area of the EBR behind a point X0 defined as the point where the fuselage starts to have shedding vortex


% Fuselage Lift
C_LB=(alpha-alpha_0B)/S *( (K*pi*D^2)/2  + eta*C_DC*(alpha-alpha_0B)*S_P_X0);


%%% WING BODY LIFT %%%
% k_wb & k_bw = factors of interference of wing-bbody for lift coefficient obtained from chart
% S_e = wing area not covered by fuselage

C_LWB = C_LB + (k_wb - k_bw)*alhpa_w*(alpha-alpha_0W)*(S_e/S); 

%%% HORIZONTAL TAIL LIFT %%% 
% a_1 is the lift slope of the tail
C_LT= a_1*(alpha*(1-d_epsilon/d_alpha) -(epsilon_0 +i_t))*(S_T/S);


%%% TOTAL LIFT %%% 
C_L=C_LW+C_LT;






end