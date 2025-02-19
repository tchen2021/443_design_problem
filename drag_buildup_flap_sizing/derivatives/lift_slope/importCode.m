% Roman Niemiec Group 5 AERO senior design Drag build up code
clc;clear;

% Import data file
DragBuildUp = importfileTEST('DragBuildUp.xlsx');
% Take individual variables from excel sheet
for i=1:height(DragBuildUp)
    assignin('base',DragBuildUp{i,2},DragBuildUp{i,3});
end

%% WING SECTION %%
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
Cl_alpha_Wing= ( (2*pi)/(2+sqrt( ((A^2 * beta_Wing^2)/2)*(1+ (tand(Lambda_half)^2)/beta_Wing^2 )+4))) /A % Function mentioned on graph 

%% TAIL SECTION %%
% Solve Graph B11A for K value
[K_Vert] = B11A(h,V,C_Vert_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Vert] = B11B(t_Vert_root,t_Vert_tip,C_Vert_root,C_Vert_tip);

%%% 5. Calculate beta
[M_Vert] = machCalc(V,h);
beta_Vert=sqrt(1-M_Vert^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Vert=(1.05/beta_Vert)*K_Vert*Cl_alpha_theory_Vert;

%%% 6. Calculate k
k_Vert=(beta_Vert*Cl_alpha_Vert)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);
% [CL_alpha_Wing] = B12(xaxis)/A;
Cl_alpha_Vert= ( (2*pi)/(2+sqrt( ((A^2 * beta_Vert^2)/2)*(1+ (tand(Lambda_Tail_half)^2)/beta_Vert^2 )+4))) /A % Function mentioned on graph 

