function [percent_power_req, percent_power_total] = engine_model(h,V,CL,CD,S)
%main engine model script
%ASSUME MAX RPM IS 2000
%   INPUTS:
        %h: altitude [ft]
        %V: speed in knots
        %S: wing area [ft^2]
%   OUTPUTS:
        %percent power required
%   STEPS:

    %calculated power required
    %check engine + prop model is able to meet that power requirement
    %calculate % power required if yes
    P_max = 1600 * 550; % lbf*ft/s
    rpm_max = 2000;
%calculate power required
    [~, ~, rho, ~] = atmosphere(h);
    V_fts = V * 1.68781; %kts to ft/s
    q = 0.5 * rho * V_fts^2;
    P_req = CD * q * S * V_fts

    %check if engine model can meet power requirement
    torque_percent = interpolateTorque(h, 0);
    [eta, ~] = interpolate_eta(V,rpm_max,h)
    
    T_max = P_max / (rpm_max * (pi/30)); %lbf-ft
    T_avi = (torque_percent/100) * T_max;

    P_avi = T_avi * (rpm_max * (pi/30)) * eta %lbf*ft/s

    if P_avi<P_req
        error('Engine model unable to provide adequate power');
    end

    %return output
    percent_power_req = (P_req/P_avi) * 100;
    percent_power_total = (P_req/P_max) * 100;

    %% plotting
    x = [-43.4 -6.3 29.7]; % temp [C]
    y = [58.5 100 100];    % percent torque
    z = [31000 12500 0];   % altitude [ft]
    n = 300;

    x_interpolated = [linspace(x(1),x(2),n), linspace(x(2),x(3),n)];
    y_interpolated = [linspace(y(1),y(2),n), linspace(y(2),y(3),n)];
    z_interpolated = [linspace(z(1),z(2),n), linspace(z(2),z(3),n)];

    T_avi = (y_interpolated./100) * T_max;
    P_avi = (((T_avi .* (rpm_max * (pi/30)) .* eta)./550) / 1600) * 100; %hp

    figure

    blue = '#2E5F7F';
    plot(P_avi, z_interpolated./1e3, 'LineWidth', 3, 'Color',blue);
    hold on
    yline(12.5, '--', 'LineWidth', 2.5)
    yline(10, '--', 'LineWidth', 2.5);
    yline(25, '--', 'LineWidth', 2.5);
    xlabel("% Cruise Power Available");ylabel("Altitude [ft x 1000]");grid on;
    fontSize_axes = 26;
    fontSize_text = 28;
    fontSize_subtitles = 28;
offset = 0.03; %offset from the horizontal line
yLineWidth = 3;

xlim([45 100]);           % X-axis limit from 0 to 0.2
ylim([0 31]);           % Y-axis limit from 0 to 1.4

% Set major ticks for grid lines
xticks(45:5:100);      % Major ticks every 0.02 on X-axis
yticks(0:5:31);       % Major ticks every 0.2 on Y-axis

% Enable grid only at major ticks
grid on;                 % Enable grid lines at major ticks
ax = gca;                % Get current axes
ax.GridLineStyle = '-';  % Solid line for major grid

% Set minor ticks without grid lines
ax.XMinorTick = 'on';           % Enable minor ticks on X-axis
ax.YMinorTick = 'on';           % Enable minor ticks on Y-axis
ax.MinorGridLineStyle = 'none'; % Turn off minor grid lines

% Define minor tick intervals
ax.XAxis.MinorTickValues = 45:1:100;  % Minor ticks every 0.005 on X-axis
ax.YAxis.MinorTickValues = 0:1:31;   % Minor ticks every 0.05 on Y-axis


% Set font size and other formatting adjustments to match style
ax.FontSize = fontSize_axes;        % Adjust font size
ax.LineWidth = 1;        % Set axis line width

end

