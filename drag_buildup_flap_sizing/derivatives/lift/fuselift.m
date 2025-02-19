function C_LB = fuselift(lf_D, alpha, alpha_0B, S, D, S_P_x0, M, l_A, l_f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% interpolate figure to get K
data_K = importfilefuseliftK("fuselift_K.csv");
 % Check that lambda is within the valid range
    if lf_D < 2 || lf_D >18
        error('the fuse length - diameter ratio must be between 2 and 18.');
    end
    
    % Extract columns
    x_values_K = data_K(:, 1);
    y_values_K = data_K(:, 2);
    
    % Interpolate y values at x_query for K
    K = interp1(x_values_K, y_values_K, lf_D, 'linear', 'extrap');

    
    % Compute the interpolated y value based on lf_D
%% interpolate figure to get eta
   data_eta = importfilefuselifteta("fuselift_eta.csv");
   %extract columns
    x_values_eta = data_eta(:, 1);
    y_values_eta = data_eta(:, 2);
    
    % Interpolate y values at x_query for K
    eta = interp1(x_values_eta, y_values_eta, lf_D, 'linear', 'extrap');

%% interpolate figure to get cdc
   data_cdc = importfilefuseliftcdc("fuselift_cdc.csv");
   %extract columns
    x_values_cdc = data_cdc(:, 1);
    y_values_cdc = data_cdc(:, 2);
    
    % Interpolate y values at x_query for K
    M_c = M*sind(alpha-alpha_0B);
    cdc = interp1(x_values_cdc, y_values_cdc, M_c , 'linear', 'extrap');

    
%% interpolate figure to get x0_lf
   data_x0_lf = importfilefuseliftx0_lf("fuselift_x0_lf.csv");
   %extract columns
    x_values_x0_lf = data_x0_lf(:, 1);
    y_values_x0_lf = data_x0_lf(:, 2);
    lA_lf = l_A/l_f;
    % Interpolate y values at x_query for K
    x0_lf = interp1(x_values_x0_lf, y_values_x0_lf, lA_lf, 'linear', 'extrap');

    
%%    

C_LB = (deg2rad(alpha-alpha_0B) ./ S) .* ( ((K.*pi.*D.^2) / (2)) + eta.*cdc .*(deg2rad(alpha)-deg2rad(alpha_0B)).*S_P_x0 ); %lift of a fuse
end

