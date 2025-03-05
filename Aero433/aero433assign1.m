%aero 433
% Global variable to store CL/CD (L/D) ratio
global CL_CD_value  % Declare global variable

%% Problem 1
options = optimoptions('fsolve', 'Display', 'off');
x0 = [3000,2000];

fun = @prob1;
x = fsolve(fun,x0,options);
% Wto1 = x(1) + (.25 * x(2));

fprintf('The maximum take-off weight for Prob 1 is %.2f lbs.\n', x(1));
fprintf('The Fuel weight for Prob 1 is %.2f lbs.\n', x(2));
fprintf('The L/D ratio for Prob 1 is %.2f\n', 6);  % L/D given as 6 for Problem 1

%% Problem 2
global CL_CD_value  % Declare the global variable in main code

fun2 = @prob2;
x1 = [3000,2000];
x2 = fsolve(fun2,x1,options);

% Wto2 = x2(1) + (.25 * x2(2));

fprintf('The maximum take-off weight for Prob 2 is %.2f lbs.\n', x2(1));
fprintf('The Fuel weight for Prob 2 is %.2f lbs.\n', x2(2));
fprintf('The L/D ratio for Prob 2 is %.2f\n', CL_CD_value);  % Print L/D value stored globally

%% Problem 3
fun3 = @prob3;
x1 = [5000,4000];  % Initial guess
x3 = fsolve(fun3,x1,options);
% Wto3 = x3(1) + (.25 * x3(2));
fprintf('The maximum take-off weight for Prob 3 is %.2f lbs.\n', x3(1));
fprintf('The Fuel weight for Prob 3 is %.2f lbs.\n', x3(2));

fprintf('The L/D ratio for Prob 3 is %.2f\n', CL_CD_value);  % Print L/D value stored globally

%% Problem 4
fun4 = @prob4;
x1 = [1000,200];  % Initial guess
x4 = fsolve(fun4,x1,options);
% Wto4 = x4(1) + (.25 * x4(2));

fprintf('The maximum take-off weight for Prob 4 is %.2f lbs.\n', x4(1));
fprintf('The Fuel weight for Prob 4 is %.2f lbs.\n', x4(2));

fprintf('The L/D ratio for Prob 4 is %.2f\n', CL_CD_value);  % Print L/D value stored globally

%% Functions

function G = prob2(y)
    global CL_CD_value  % Declare global variable inside the function
    Np = 0.85;
    SFC = 0.4;
    Range = 1050 / 325;
    Wpay = (4 * 240) + 80;
    S = 13.5 * 10.764;
    rho = 0.0023769;  % slug/ft^3
    Cd0 = 0.025;
    e = 0.85;
    AR = 10.1;
    Vcruise = 140 * 1.688;  % ft/s
    
    % Lift coefficient and drag coefficient
    CL = (2 * y(1)) / (rho * Vcruise^2 * S);
    Cd = Cd0 + (CL^2) / (pi * e * AR);
    CL_CD_value = CL / Cd;  % Store CL/CD ratio globally
    
    G(1) = (((Np * CL_CD_value) / SFC) * log((y(1)) / (y(1) - y(2)))) - Range;
    G(2) = (0.5 * y(1)) - Wpay - y(2);
end

function h = prob3(y)
    global CL_CD_value  % Declare global variable inside the function
    Np = 0.85;
    SFC = 0.4;
    Range = 1050 / 325;
    Wpay = (4 * 240) + 80;
    S = 13.5 * 10.764;
    rho = 0.0023769;  % slug/ft^3
    Cd0 = 0.025;
    e = 0.85;
    AR = 10.1;
    Vcruise = 180 * 1.688;  % ft/s
    
    % Lift coefficient and drag coefficient
    CL = (2 * y(1)) / (rho * Vcruise^2 * S);
    Cd = Cd0 + (CL^2) / (pi * e * AR);
    CL_CD_value = CL / Cd;  % Store CL/CD ratio globally
    
    h(1) = (((Np * CL_CD_value) / SFC) * log((y(1)) / (y(1) - y(2)))) - Range;
    h(2) = (0.5 * y(1)) - Wpay - y(2);
end

function i = prob4(y)
    global CL_CD_value  % Declare global variable inside the function
    Np = 0.85;
    SFC = 0.4;
    Range = 1050 / 325;
    Wpay = (6 * 240) + 100;
    S = 13.5 * 10.764;
    rho = 0.0023769;  % slug/ft^3
    Cd0 = 0.025;
    e = 0.85;
    AR = 10.1;
    Vcruise = 180 * 1.688;  % ft/s
    
    % Lift coefficient and drag coefficient
    CL = (2 * y(1)) / (rho * Vcruise^2 * S);
    Cd = Cd0 + (CL^2) / (pi * e * AR);
    CL_CD_value = CL / Cd;  % Store CL/CD ratio globally
    
    i(1) = (((Np * CL_CD_value) / SFC) * log((y(1)) / (y(1) - y(2)))) - Range;
    i(2) = (0.5 * y(1)) - Wpay - y(2);
end

function F = prob1(x)
    Np = 0.85;
    L_D = 6;  % Given L/D
    SFC = 0.4;
    Range = 1050 / 325;
    Wpay = (4 * 240) + 80;
    
    F(1) = (((Np * L_D) / SFC) * log((x(1)) / (x(1) - x(2)))) - Range;
    F(2) = (0.5 * x(1)) - Wpay - x(2);
end


