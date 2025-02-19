function [alpha_0W] = zeroliftalpha(iprime_r,epsilon,alpha_0W_root, AR, lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%interpolate figure to get J
data = importfileZeroLiftAlpha("zero-lift-alpha.csv");
 % Check that lambda is within the valid range
    if lambda < 0.4 || lambda > 0.6
        error('Taper ratio must be between 0.4 and 0.6.');
    end
    
    % Extract columns
    x_values = data(:, 1);
    y_0_4 = data(:, 2);
    y_0_6 = data(:, 3);
    
    % Interpolate y values at x_query for lambda
    y_0_4_interp = interp1(x_values, y_0_4, AR, 'linear', 'extrap');
    y_0_6_interp = interp1(x_values, y_0_6, AR, 'linear', 'extrap');
    
    % Compute the interpolated y value based on lambda
    J = y_0_4_interp + (lambda - 0.4) / (0.6 - 0.4) * (y_0_6_interp - y_0_4_interp);

alpha_0W = alpha_0W_root - iprime_r + J*epsilon;
end

