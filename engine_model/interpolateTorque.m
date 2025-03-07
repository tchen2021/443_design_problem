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

%% plotting 
blue = '#2E5F7F';
plot(y_interpolated, z_interpolated./1e3, 'LineWidth', 3, 'Color',blue);
hold on
yline(12.5, '--', 'LineWidth', 2.5)
xlabel("% Torque");ylabel("Altitude [ft x 1000]");grid on;
fontSize_axes = 26;
fontSize_text = 28;
fontSize_subtitles = 28;
offset = 0.03; %offset from the horizontal line
yLineWidth = 3;

xlim([55 100]);           % X-axis limit from 0 to 0.2
ylim([0 35]);           % Y-axis limit from 0 to 1.4

% Set major ticks for grid lines
xticks(55:5:100);      % Major ticks every 0.02 on X-axis
yticks(0:5:35);       % Major ticks every 0.2 on Y-axis

% Enable grid only at major ticks
grid on;                 % Enable grid lines at major ticks
ax = gca;                % Get current axes
ax.GridLineStyle = '-';  % Solid line for major grid

% Set minor ticks without grid lines
ax.XMinorTick = 'on';           % Enable minor ticks on X-axis
ax.YMinorTick = 'on';           % Enable minor ticks on Y-axis
ax.MinorGridLineStyle = 'none'; % Turn off minor grid lines

% Define minor tick intervals
ax.XAxis.MinorTickValues = 55:1:100;  % Minor ticks every 0.005 on X-axis
ax.YAxis.MinorTickValues = 0:1:35;   % Minor ticks every 0.05 on Y-axis

% Label axes
%xlabel('C_D');           % Label for X-axis
%ylabel('C_L', Rotation= 0);           % Label for Y-axis

% Set font size and other formatting adjustments to match style
ax.FontSize = fontSize_axes;        % Adjust font size
ax.LineWidth = 1;        % Set axis line width


end

