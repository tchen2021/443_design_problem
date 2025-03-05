
clear; clc; close all;
format compact;

%% L/D ratios and SFC (What we are changing)
LD_coarse = linspace(10, 20, 4);  % Lift-to-drag ratios
SFC_coarse = linspace(0.3, 0.6, 5);  % Specific Fuel Consumption values

% Refine the grid using a custom refinement function
nref = 3;%%?
LD = refvec(LD_coarse, nref);
SFC = refvec(SFC_coarse, nref);

EmptyRatio = .52;
%empty rate ratio in book .54

% Preallocation
solution2 = zeros(length(LD), length(SFC), 2); 
Mff = zeros(9,13);
W_f = zeros(9,13);

% Set options for fsolve to suppress output
options = optimoptions('fsolve', 'Display', 'none');

% Nested for loop to solve for each LD and SFC

for i = 1:length(LD)

    for j = 1:length(SFC)

        payload = (150 * (175 + 30)) + (5 * (175 + 30));  % Payload calculation (including crew weight)
        E = 1;
        R = (1500);  % Range in nautical miles converted (flight credit considered)
        V = 473;  % Velocity in knots (cruise from phase 5)

        % Weight ratios for Transport Jets
        W1_T0 = .990; % Start up
        W2_1 = .990; % Taxi
        W3_2 = .995; % Take-off
        W4_3 = .980; % Climb

        W7_6 = .990; % Descent 
        
        W9_8 = .992; % Shutdown

        % Define the system of equations with LD(i) and SFC(j)
        % 1 = W4/5 2 = W5/6 3 = W7/8
        F = @(W) [
            
            (LD(i) * V / SFC(j)) * log(W(1)^-1) - R;  % Range equation for cruise (Phase 5)
            (18 / .6) * log(W(2)^-1) - E;             % Endurance equation for loiter (Phase 6)
            % L/D = 18, SFC = .6
            (10 * 250 / .9) * log(W(3)^-1) - 100;     % Range equation for alternate (Phase 8)
            % L/D = 10, SFC = .6, V = 250 knots
            ];

        % Initial guess for W_to and W_f
        initial_guess = [.1, .11, .12];  % We are working with small fraction ratios
        % Because the functions work with the natural logs, initial guesses
        % CANNOT BE 0 OR ALL THE SAME

        % Solve the system using fsolve
        solution = fsolve(F, initial_guess, options);

        Mff(i,j) = (W9_8 * solution(3) * W7_6 * solution(2) * solution(1) * W4_3 * W3_2 * W2_1 * W1_T0);
        W_f(i,j) = (1 - Mff(i,j));  

        %  % Display the result
        % disp(['Solution for LD = ', num2str(LD(i)), ', SFC = ', num2str(SFC(j))]);
        % disp(['W_f ratio to Wto = ', num2str(W_f(i,j))]);

        %% Weight Ratio (Predetermined ratio)
        G = @(Wto) ((Wto - payload - (W_f(i,j) * Wto) - (Wto * .005)) / Wto) - EmptyRatio;
        
        initial_guess2 = 110000;  % 110,000 lbs is the initial guess in the book

        % Solve the system using fsolve with optionsW
        solution2(i,j) = fzero(G, initial_guess2);

        

    end
end

% Define offset for the carpet plot
offset = 1;

% Create the carpet plot
figure(1)
carpet(LD, SFC, solution2(:, :, 1)', offset, nref, 'b', 'r');

% Add labels for vertical lines (corresponding to LD axis)
carpetlabel(LD, SFC, solution2(:, :, 1)', offset, nref, .2, 0, 0.1, 0.1);

% Add labels for horizontal lines (corresponding to SFC axis)
carpetlabel(LD, SFC, solution2(:, :, 1)', offset, nref, 0, -1, 0, 0.5);

%%
% Add labels for each L/D line at the top of the plot (highest SFC value)
% Loop through each L/D value 
for i = 1:length(LD)
    % Check if the current L/D value is part of the displayed data
    if ismember(LD(i), LD_coarse)  % Check if LD(i) is in the coarse values
        % Label at the highest SFC value (top of the plot)
        carpettext(LD, SFC, solution2(:, :, 1)', offset, LD(i), SFC(end), ...
            ['L/D = ', num2str(LD(i))], 0.05, 0.05, 'FontSize', 8, 'Color', 'g');
    end
end

% Add labels for each SFC line on the right side (highest L/D value)
for j = 1:length(SFC)
    % Check if the current SFC value is part of the displayed data
    if ismember(SFC(j), SFC_coarse)  % Check if SFC(j) is in the coarse values
        % Label at the highest L/D value (right side of the plot)
        carpettext(LD, SFC, solution2(:, :, 1)', offset, LD(end), SFC(j), ...
            ['SFC = ', num2str(SFC(j))], 0.05, 0.05, 'FontSize', 8, 'Color', 'b');
    end
end


