% Create a new figure
figure;

% Define axes limits
xlim([0 0.2]);           % X-axis limit from 0 to 0.2
ylim([0 1.4]);           % Y-axis limit from 0 to 1.4

% Set major ticks for grid lines
xticks(0:0.02:0.2);      % Major ticks every 0.02 on X-axis
yticks(0:0.2:1.4);       % Major ticks every 0.2 on Y-axis

% Enable grid only at major ticks
grid on;                 % Enable grid lines at major ticks
ax = gca;                % Get current axes
ax.GridLineStyle = '-';  % Solid line for major grid

% Set minor ticks without grid lines
ax.XMinorTick = 'on';           % Enable minor ticks on X-axis
ax.YMinorTick = 'on';           % Enable minor ticks on Y-axis
ax.MinorGridLineStyle = 'none'; % Turn off minor grid lines

% Define minor tick intervals
ax.XAxis.MinorTickValues = 0:0.005:0.2;  % Minor ticks every 0.005 on X-axis
ax.YAxis.MinorTickValues = 0:0.05:1.4;   % Minor ticks every 0.05 on Y-axis

% Label axes
xlabel('C_D');           % Label for X-axis
ylabel('C_L');           % Label for Y-axis

% Set font size and other formatting adjustments to match style
ax.FontSize = 10;        % Adjust font size
ax.LineWidth = 1;        % Set axis line width

% Add title or any other annotations if desired
title('Adjusted Figure Setup with Major Grid Lines and Minor Ticks');
