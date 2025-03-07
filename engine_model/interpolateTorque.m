function torque_percent = interpolateTorque(altitude, deviation)
%interpolates engine data based on PT6-68B maximum cruise performance p343
%of pdf of AT-6 POH

%%%%%%%%%INPUTS
%altitude [ft]
%velocity [kts]
%S        [ft^2]

%%%%%%%%%%OUTPUTS
%torque_percent: percent torque

%define data points from engine model
x = [-43.4 -6.3 29.7]; % temp [C]
y = [58.5 100 100];    % percent torque
z = [31000 12500 0];   % altitude [ft]
n = 300;

% Call the atmosphere function to get atmospheric properties
[~, T_rankine, ~, ~] = atmosphere(altitude);

% Convert temperature from Rankine to Celsius
T_celsius = ((T_rankine - 491.67) * 5/9) + deviation;

% if T_celsius<x(1) || T_celsius>x(3)
%     error('temp from alt and deviation outside of range')
% end

if altitude>z(1) || altitude<z(3)
    error('alt outside of range')
end

%interpolate defined data points

x_interpolated = [linspace(x(1),x(2),n), linspace(x(2),x(3),n)];
y_interpolated = [linspace(y(1),y(2),n), linspace(y(2),y(3),n)];
z_interpolated = [linspace(z(1),z(2),n), linspace(z(2),z(3),n)];

%interpolate by stdatm temperature from alt input
% margin = (abs(x(1))-abs(x(2))) / n;
% 
% valid_indices = find(abs(x_interpolated - T_celsius) <= margin);
% 
% % Find the index of the closest value within the valid indices
% [~, min_idx] = min(abs(x_interpolated(valid_indices) - T_celsius));
% 
% % Return the corresponding index in the original array
% idx_interpolation = valid_indices(min_idx);

%interpolate directly by altitude
margin = (abs(z(1))-abs(z(2))) / n;

valid_indices = find(abs(z_interpolated - altitude) <= margin);

% Find the index of the closest value within the valid indices
[~, min_idx] = min(abs(z_interpolated(valid_indices) - altitude));

% Return the corresponding index in the original array
idx_interpolation = valid_indices(min_idx);

torque_percent = y_interpolated(idx_interpolation);
end

