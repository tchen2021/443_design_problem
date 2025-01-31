function [K] = B11A(h,V,C_root)
%Code figures out linear interpolation of graph B11A
% Starting with sutherlands then solving for Re, the code uses that value
% and the linear interpolation between Re #s from 1,000,000 - 10,000,000 Re
% to figure out what the value of K is.

%% 1. Sutherlands Equation 
% Get variables
Sutherlands_air_const=1.716E-5; 
T_0=518.7; % [R]
[p, T, rho, a] = atmosphere(h); % get temp in [R]
Sutherlands= 3.62E-7 * ((T/T_0)^1.5)* ((T_0 + 198.72)/(T+198.72));
%Sutherland=	0.000012276873535233*((.555*T_0 +120)/(.555*T + 120))*((T/T_0)^1.5)

%% Reynolds number
Re= (rho*V*C_root)/Sutherlands;

%% Linear Interpolate from graphs
% load digitized graph data
B111AE7 = importB11A("B11-1A_E7.csv");
B111AE6 = importB111AE6("B11-1A_E6.csv");
halfPhiPrime=10.125;

if Re> 1E6 && Re<1E7
    [~,Index106] = min(abs(tand(halfPhiPrime)-B111AE6(:,1))); % Find location along 10E6 line from LE Angle
    [~,Index107] = min(abs(tand(halfPhiPrime)-B111AE7(:,1))); % Find location along 10E7 line from LE Angle
    K=-(((1E7-Re)/(1E7-1E6))*(B111AE7(Index107,2)-B111AE6(Index106,2)) - B111AE7(Index107,2)); % Linear interpolation
end

% Solve for other inputs

end % end of function