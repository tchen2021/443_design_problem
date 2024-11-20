% Lift and Drag coefficients
CLcr = [1.98E-01, 1.24E-01, 2.64E-01, 3.02E-01, 2.52E-01, 2.88E-01, 1.50E-01, 1.74E-01];
CDcr = [2.79E-02, 2.53E-02, 3.56E-02, 3.60E-02, 3.83E-02, 4.70E-02, 2.07E-02, 1.69E-02];

% Unique names for each data point (excluding P5, P11, and P13)
point_names = {'P1', 'P2', 'P3', 'P4', 'P6', 'Bronco', 'P10', 'P12'};

% Updated indices for jets and propellers
jet_indices = [2, 4];
prop_indices = setdiff(1:length(CLcr), jet_indices);

% OV-10 Bronco optimized drag polar parameters
C_D0 = 0.023;     % Drag coefficient at zero lift
e = 0.85;         % Oswald efficiency factor
AR = 6;         % Aspect ratio

% Bronco optimized drag polar
C_L_bronco = linspace(0, 1.1, 40);
C_D_bronco = C_D0 + (C_L_bronco.^2) / (pi * e * AR);

% This is for the bad case
C_D0_bad = 0.038;    % Higher drag coefficient at zero lift
e_bad = 0.85;        % Lower Oswald efficiency factor
AR_bad = 4.5;        % Lower aspect ratio

% Bronco bad case drag polar
C_D_bronco_bad = C_D0_bad + (C_L_bronco.^2) / (pi * e_bad * AR_bad);

% Define the line for L/D = 
target_L_D = 9;
C_L_tangent = linspace(0, 1.1, 40);
C_D_tangent = C_L_tangent / target_L_D;

% Define another line for L/D = 
target_L_D1 = 13;
C_L_tangent1 = linspace(0, 1.1, 40);
C_D_tangent1 = C_L_tangent1 / target_L_D1;

CD01 = .023;
CD02 = .038;

AR1 = 6.609;
AR2 = 4;
e1 = .85;
k = 1/(pi*AR1*e1);
k1 = 1/(pi*AR2*e1);
Line1 = sqrt(CD01/k);
Line2 = sqrt(CD02/k1);


figure;
hold on;

% Plot and label propeller points
for i = prop_indices
    plot(CDcr(i), CLcr(i), 'bo'); % Plot propeller points in blue
    text(CDcr(i), CLcr(i), point_names{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Plot and label jet points
for i = jet_indices
    plot(CDcr(i), CLcr(i), 'ro'); % Plot jet points in red
    text(CDcr(i), CLcr(i), point_names{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Plot optimized Bronco drag polar line
plot(C_D_bronco, C_L_bronco, 'k-', 'LineWidth', 1.5, 'DisplayName', 'OV-10 Bronco Optimized');

% Plot bad case Bronco drag polar line
plot(C_D_bronco_bad, C_L_bronco, 'k-', 'LineWidth', 1.5, 'DisplayName', 'OV-10 Bronco Bad Case');

% Plot the L/D = 9.6 tangent line
plot(C_D_tangent, C_L_tangent, 'g--', 'LineWidth', 1.5, 'DisplayName', 'L/D = 9.6');

% Plot the L/D = 11.5 tangent line
plot(C_D_tangent1, C_L_tangent1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'L/D = 11.5');
%scatter(CD_dat, CL_dat, 'filled', MarkerFaceColor = "#D95319", MarkerEdgeColor = "#D95319")
%plot(CDtot_bad, CL, 'k')
% Plot settings
yline(Line1)
yline(Line2,'g')
xlabel('Drag Coefficient (C_D)');
ylabel('Lift Coefficient (C_L)');
title('Optimized vs. Bad Case Drag Polar for OV-10 Bronco');
hold off;