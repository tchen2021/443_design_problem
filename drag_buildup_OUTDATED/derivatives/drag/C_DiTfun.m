function [C_DiT] = C_DiTfun(C_LT, A_T, S, S_T, t_c, beta, cf_C)
%
    % Import data sets
    data1 = importfiledrag_tail_delta1("drag_tail_delta1.csv"); % Call function to get first data set
    data2 = importfiledrag_tail_delta2("drag_tail_delta2.csv"); % Call function to get second data set

    % Extract x and y data from the first data set
    % Assume data1 has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data1(:, 1); % First column is x data
    y1_1 = data1(:, 2); % Second column is y data for curve 1
    y1_2 = data1(:, 3); % Third column is y data for curve 2
    y1_3 = data1(:, 4); % Fourth column is y data for curve 3

    % Extract x and y data from the second data set
    % Assume data2 has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x2 = data2(:, 1); % First column is x data
    y2_1 = data2(:, 2); % Second column is y data for curve 1
    y2_2 = data2(:, 3); % Third column is y data for curve 2
    y2_3 = data2(:, 4); % Fourth column is y data for curve 3

    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, cf_C, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, cf_C, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, cf_C, 'linear', 'extrap');

    % Linearly interpolate between the three curves based on t_c

    if t_c >=0.12 && t_c<=0.21
        % Interpolate between curve 1 and curve 2
        tail_delta1 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
    elseif t_c<0.12 || t_c >0.3
        error("t/c must be between 0.12 and 0.3");
    else
        % Interpolate between curve 2 and curve 3
        tail_delta1 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
    end

    % Perform interpolation for the second data set using A_T as the query point
    y2_1_interp = interp1(x2, y2_1, beta, 'linear', 'extrap');
    y2_2_interp = interp1(x2, y2_2, beta, 'linear', 'extrap');
    y2_3_interp = interp1(x2, y2_3, beta, 'linear', 'extrap');

    % Linearly interpolate between the three curves based on t_c
    if t_c <= 0.3 && t_c>=0.21
        % Interpolate between curve 1 and curve 2
        tail_delta2 = y2_1_interp + (y2_2_interp - y2_1_interp) * t_c / 0.5;
    elseif t_c<0.12 || t_c > 0.3
        error("t/c must be between 0.12 and 0.3")
    else
        % Interpolate between curve 2 and curve 3
        tail_delta2 = y2_2_interp + (y2_3_interp - y2_2_interp) * (t_c - 0.5) / 0.5;
    end

    % Output the interpolated results
    % You can manually define the final output using tail_delta1 and tail_delta2
    C_DiT = ((C_LT^2) / (pi * A_T)) * (1+tail_delta1+tail_delta2) * (S_T/S); %induced drag of the horizontal tail
end