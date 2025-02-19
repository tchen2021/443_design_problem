function [C_L_alpha] = C_L_alpha(t_root,t_tail,C_root,C_tip,h,V, M,AR,Lambda_half)
%Calculates C_L_alpha using interpolation algorithms 
%% interpolating for K
K = B11A(h, V, C_root);


%% Interpolating for  Cl_alpha theory
Cl_alpha_theory = B11B(t_root,t_tail,C_root,C_tip);

%% calculating Prandtl-Galuert compressibility factor
beta = sqrt(1-M^2);

%% calcuating lift slope of 2D wing alone
C_l_alpha = (1.05/beta) * K * Cl_alpha_theory;

%% calculating 3D lift slope
C_L_alpha_A = B121(AR, K, beta, Lambda_half);

C_L_alpha = C_L_alpha_A * AR;


end % end of function