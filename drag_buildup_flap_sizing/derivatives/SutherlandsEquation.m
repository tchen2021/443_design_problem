function [Sutherlands] = SutherlandsEquation(h)
%This function solves for sutherlands value (mu) to then be used for Re for
%air
% Use sutherlands constant for air
Sutherlands_air_const=1.716E-5; 
T_0=518.7; % [R]
[p, T, rho, a] = atmosphere(h); % get temp in [R]
Sutherlands= 3.62E-7 * ((T/T_0)^1.5)* ((T_0 + 198.72)/(T+198.72));
end