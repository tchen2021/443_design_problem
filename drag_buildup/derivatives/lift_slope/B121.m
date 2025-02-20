function [C_L_a_A] = B121(AR, K, beta, Lambda_half)
%Code figures out linear interpolation of graph B11A
% Starting with sutherlands then solving for Re, the code uses that value
% and the linear interpolation between Re #s from 1,000,000 - 10,000,000 Re
% to figure out what the value of K is.

%% 
%obtaining x value
x = (AR/K) * (beta^2 + tand(Lambda_half)^2)^0.5;
%% Linear Interpolate from graphs
% load digitized graph data
if x> 16 && x < 0
    error("x value needs to be between 0 and 16") 
end
data = importfileB121("B12-1");


% Extract columns
x_values = data(:, 1);
y_values = data(:, 2);

% Interpolate y values at x_query for K
C_L_a_A = interp1(x_values, y_values, x, 'linear', 'extrap');


    




end % end of function