% Roman Niemiec Group 5 AERO senior design Drag build up code
clc;clear;

% Import data file
DragBuildUp = importfileTEST('DragBuildUp.xlsx');
% Take individual variables from excel sheet
for i=1:height(DragBuildUp)
    assignin('base',DragBuildUp{i,2},DragBuildUp{i,3});
end

% Solve Graph B11A for K value
[K] = B11A(h,V,C_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory] = B11B(t_root,t_tail,C_root,C_tip);

%%% 4. Calculate Cl_alpha
% Initialize variables
beta
Cl_alpha=(1.05/beta)*K*Cl_alpha_theory