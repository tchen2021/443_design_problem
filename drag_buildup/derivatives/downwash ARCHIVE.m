function [deltaepsilon_deltaalpha] = downwash(Lambda,AR,lambda, tailfactor, armfactor,h_H,b,l_H)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
%interpolate figure to get K_A
data_K_A = importfileB51("B51.csv");
 % Check that lambda is within the valid range
    if AR >= 10  || AR <= 0
        error('Aspect ratio must be between 0 and 10.');
    end
    
    % Extract columns
    x_K_A = data_K_A(:, 1);
    y_K_A = data_K_A(:, 2);
    
    
    % Interpolate y values at x_query for lambda
    K_A_interp = interp1(x_K_A, y_K_A, AR, 'linear', 'extrap');
    %y_K_A_interp = interp1(x_K_A, y_K_A, AR, 'linear', 'extrap');
    
    % Compute the interpolated y value based on lambda
    %J = x_K_A_interp + (lambda - 0.4) / (0.6 - 0.4) * (y_K_A_interp - x_K_A_interp);

%% interpolate figure to get K_B
%interpolate figure to get K_A
data_K_lambda = importfileB52("B52.csv");
 % Check that lambda is within the valid range
    if lambda >=1  || lambda <= 0
        error('Taper ratio must be between 0 and 1.');
    end
    
    % Extract columns
    x_K_lambda = data_K_lambda(:, 1);
    y_K_lambda = data_K_lambda(:, 2);
    
    
    % Interpolate y values at x_query for lambda
    K_lambda_interp = interp1(x_K_lambda, y_K_lambda, lambda, 'linear', 'extrap');

%%
%interpolate figure to get K_H
data_K_H = importfileB53("B53.csv");
 % Check that tail factor abs(2h_H/b) is within the valid range
    if tailfactor >=1  || tailfactor <= 0
        error('Tail Factor must be between 0 and 1.');
    end
    
    % Extract columns
    x_K_tailfactor = data_K_H(:, 1);
    y_K_tailfactor_0_4 = data_K_H(:, 2);
    y_K_tailfactor_0_6 = data_K_H(:, 3);
    y_K_tailfactor_0_8 = data_K_H(:, 4);
    
    % Interpolate y values at x_query for lambda
    y_K_0_4_interp = interp1(x_K_tailfactor, y_K_tailfactor_0_4, AR, 'linear', 'extrap');
    y_K_0_6_interp = interp1(x_K_tailfactor, y_K_tailfactor_0_6, AR, 'linear', 'extrap');
    y_K_0_8_interp = interp1(x_K_tailfactor, y_K_tailfactor_0_8, AR, 'linear', 'extrap');
    
    
    if armfactor <=0.6  || armfactor >= 0.4

    % Compute the interpolated y value based on lambda
    K_H = y_K_0_4_interp + (armfactor - 0.4) / (0.6 - 0.4) * (y_K_0_6_interp - y_K_0_4_interp);
    
    elseif armfactor <=0.8  || armfactor >= 0.6

    % Compute the interpolated y value based on lambda
    K_H = y_K_0_6_interp + (armfactor - 0.6) / (0.8 - 0.6) * (y_K_0_8_interp - y_K_0_6_interp);
    
    else 
        error('Arm Factor must be between 0.4 and 0.8.');
    end


%%% USE THE EQUATION
K_H= (1-abs(h_H/b))/ (((2*l_H)/b)^1.5)

%% final calculation
deltaepsilon_deltaalpha = 4.44 * (K_A_interp * K_lambda_interp * K_H * cosd(Lambda)^0.5)^1.19;

