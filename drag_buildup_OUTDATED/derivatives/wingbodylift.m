function C_LWB = wingbodylift(C_LB, a_W, alpha, alpha_0W, Se_S, D_b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% interpolate figure to get K_W and K_B
data = importfilewingbodylift("wingbody_lift.csv");
 % Check that lambda is within the valid range
    if  D_b< 0 || lf_D >1
        error('the fuse diameter/wing span ratio must be between 0 and 1.');
    end
    
    % Extract columns
    x_values = data(:, 1);
    y_values_W = data(:, 3);
    y_values_B = data(:, 2);
    
    % Interpolate y values at x_query for K
    K_W = interp1(x_values, y_values_W, D_b, 'linear', 'extrap');
    K_B = interp1(x_values, y_values_B, D_b, 'linear', 'extrap');

    
    % Compute the interpolated y value based on lf_D

%%    

C_LWB = C_LB + (K_W-K_B).*a_W.*(alpha-alpha_0W).*(Se_S);     %total lfit of wing-body system, alpha in rad


