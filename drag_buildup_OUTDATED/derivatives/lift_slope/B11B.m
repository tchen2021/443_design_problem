function [Cl_alpha_theory] = B12(t_root,t_tail,C_root,C_tip)
%Code figures out based on x value what the y value that matches with it
%should be

%% 1. Solve for x axis
% get average values of thickness and chord
t=(t_root+t_tail)/2;
C=(C_tip+C_root)/2;
wingThicknessRatio=t/C_root;

%% Fnd point on graph (B11B)
% load digitized graph data
B111B_input = importfileB11B("B11-1B.csv");
[~,Index] = min(abs(wingThicknessRatio-B111B_input(:,1))); % Find location along 10E7 line from LE Angle
Cl_alpha_theory=B111B_input(Index,2);

end % end of function