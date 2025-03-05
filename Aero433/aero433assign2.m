%aero 433 2
%influence over SFC , L/D and empty weight factor
% will go with .4-.65 SFC range as seen in chatgpt and gudmundson book
% L_D boeing 747 takeoff 7-10, climb 10- 15, cruise 17-18
% L_D Airbus A320 takeoff 6-8, climb 10-12, cruise 15-16
%DC-9-80 takeoff 5-7, climb 9-11, cruise 14-15,
format compact;
% empty weight factor boeing 747 .5-.52, Aribus 1320 .45-.5, DC-9-80 .5-.55

%using Wpay = 40,000 just in case we need to go heavier
% Np boeing .65-.7, Airbus A320 .65-.75, DC-9-80 .55-.6, C-17, .6-.7 
clc
clear
options = optimoptions('fsolve', 'Display', 'off');
L_Dcoarse = linspace(15, 17, 5);  % Given L/D
SFCcoarse = linspace(0.4, 0.65,5); % Specific fuel consumption

nref = 3;
LD = refvec(L_Dcoarse,nref);
SFC = refvec(SFCcoarse,nref);
% Initial guess for x
x0 = [130000, 50000];

% Pre-allocate array to store solutions for each combination of L/D and SFC
x_solutions = zeros(2,length(L_Dcoarse), length(SFCcoarse));

% Loop through each value of L/D and SFC
for i = 1:length(LD)
    for j = 1:length(SFC)
         L_D = L_Dcoarse(i);
         SFC = SFCcoarse(j);

        % Call fsolve for each combination of L/D and SFC
        x_solution = fsolve(@(x) prob1(x, L_D, SFC), x0,options);

        % Store the result in x_solutions array
        x_solutions(:, i, j) = x_solution;

    end
end
% 
% 
% % Define offset for the carpet plot
% offset = 1;
% % 
% % % Create the carpet plot
%  figure(1)
%  carpet(LD, SFC, x_solutions(:, :, 1)', offset, nref, 'b', 'r');
% % 
% 





function F = prob1(x, L_D, SFC)

    V = 473*1.68; % Vel
    Range = 1600/ 325; % Nm
  Wpay = (150 * (175 + 30)) + (2 * (175 + 30));  % Payload calculation
  weightfrac = .53;
    % Solve for a single pair of L_D and SFC
    F(1) = (L_D .* V .* (1./SFC) .* (log(x(1) ./ (x(1) - x(2))))) - Range;
    F(2) = (weightfrac .* x(1)) - Wpay - x(2);

end





