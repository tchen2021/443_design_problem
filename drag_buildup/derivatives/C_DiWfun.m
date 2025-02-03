function C_DiW = C_DiWfun(lambda, Delta_quarter, C_LW, AR)
    % Import data sets
    data1 = importfiledragdelta1("drag_delta1"); % Call function to get first data set
    data2 = importfiledragdelta2("drag_delta2"); % Call function to get second data set

    x1 = data1(:, 1); % First column is x data
    y1_4 = data1(:, 2); % Second column is y data for curve 4
    y1_6 = data1(:, 3); % Third column is y data for curve 6
    y1_8 = data1(:, 4); % Fourth column is y data for curve 8

    % Extract x and y data from the second data set
    x2 = data2(:, 1); % First column is x data
    y2 = data2(:, 2); % Second column is y data

    % Perform interpolation for the first data set using lambda as the query point
    y1_4_interp = interp1(x1, y1_4, lambda, 'linear', 'extrap');
    y1_6_interp = interp1(x1, y1_6, lambda, 'linear', 'extrap');
    y1_8_interp = interp1(x1, y1_8, lambda, 'linear', 'extrap');

    % Linearly interpolate between the three curves based on AR
    if AR <= 5
        % Interpolate between curve 4 and curve 6
        delta1 = y1_4_interp + (y1_6_interp - y1_4_interp) * (AR - 4) / 2;
    else
        % Interpolate between curve 6 and curve 8
        delta1 = y1_6_interp + (y1_8_interp - y1_6_interp) * (AR - 6) / 2;
    end
    % Perform interpolation for the second data set using Delta_quarter as the query point
    delta2 = interp1(x2, y2, Delta_quarter, 'linear', 'extrap');

    % Combine the interpolated results (e.g., sum, average, or other operation)
    C_DiW = ((C_LW^2)/(pi*AR)) * (1+ delta1 + delta2); % Example: sum of the two interpolated values
end