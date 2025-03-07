function torque_percent = interpolateTorque_defunct(altitude, velocity, CL, CD, S, eta, deviation)
%%%%%%%%%INPUTS
%altitude [ft]
%velocity [kts]
%S        [ft^2]

%%%%%%%%%%OUTPUTS
%torque_percent: percent torque

%%%%%%NOTE:
%FUNCTION CURRENTLY NOT WORKING

    % Call the atmosphere function to get atmospheric properties
    [~, T_rankine, ~, ~] = atmosphere(altitude);
    
    % Convert temperature from Rankine to Celsius
    T_celsius = ((T_rankine - 491.67) * 5/9) + deviation;
    
    % Import the digitized data
    data = importfileGeneral("pt6A-68_max_cruise_power_final.csv", 18);
    
    % Extract the x (temperature) and y (torque %) data
    x_data = data(:, 1);  % Temperature data (x-axis)
    y_data_altitude = data(:, 2:10);  % Torque data for constant altitude lines
    y_data_temp_offset = data(:, 11:18);  % Torque data for constant temperature offset lines
    
    % Define the altitude and temperature offset values
    altitude_values = [0, 4000, 8000, 12000, 16000, 20000, 24000, 28000, 31000];  % in feet
    temp_offset_values = [37	30	20	10	0	-10	-20	-30];  % in Celsius
    
    % Create a 2D grid for interpolation
    % Combine the altitude and temperature offset data into a single 2D grid
    [X, Y] = meshgrid(temp_offset_values, altitude_values);
    
    % Reshape the torque data into a 2D grid
    % Each row corresponds to a temperature offset, and each column corresponds to an altitude
    torque_grid = zeros(length(altitude_values),length(temp_offset_values));

  
    for i = 1:length(altitude_values)
        for j = 1:length(temp_offset_values)
            %obtain the range of values from the ith and jth const. alt and
            %const. temp offset values
            min_alt = min(y_data_altitude(:,i));
            max_alt = max(y_data_altitude(:,i));
            min_temp_offset = min(y_data_temp_offset(:,j));
            max_temp_offset = max(y_data_temp_offset(:,j));

            % Interpolate torque for each combination of altitude and temperature offset

            torque_grid(i, j) = interp1(x_data, y_data_altitude(:, j), T_celsius, 'linear', 'extrap');
            if torque_grid(i,j)<min_alt || torque_grid(i,j)>max_alt || torque_grid(i,j)<min_temp_offset || torque_grid(i,j)>max_temp_offset %check that the interpolate value exists for the ith and jth combination of const. alt and const. temp offset
                error('combination of altitude and temp offset at query point does not exist');
            end
        end
    end

      
    
    % Perform 2D interpolation to find the torque at the given altitude and temperature deviation
    torque_percent = interp2(X, Y, torque_grid, altitude, deviation, 'linear');
    
    % Ensure the torque percentage is within the range 30 to 100
    torque_percent = max(30, min(100, torque_percent));
end