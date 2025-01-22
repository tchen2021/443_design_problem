%% Class 2 estimation
clear
clc

%% Weight Estimation Raymer Methods
% Date Edited: 01/15/25

% Authors: Michael Chen
%          Jacob Quiray

% Inputs:

% MTOW
% Wing Parameters
% Fuselage Parameters

% Outputs:

% Weight of each component
% CG of each component as a % of MAC

%% Validation Test

%% Wing Weight

% Fighter/Attack
% Define the variables Embraer
    Embraer.K_dw = 1; % Wing design factor
    Embraer.K_vs = 1; % Varible sweep factor
    Embraer.W_dg = 12500; % Design gross weight (lbs or other units as needed)
    Embraer.Nz = 3.75; % Ultimate Load factor %ASSUMED
    Embraer.S_w = 313; % Wing reference area (sq ft or other units)
    Embraer.A = 8.07; % Aspect ratio
    Embraer.t_c_root = .16; % Thickness-to-chord ratio at the root
    Embraer.lambda = 0.55; % Taper ratio
    Embraer.Lambda = 0; % Sweep angle (in radians)
    %Embraer.S_csw = 2*11.55; % Control surface area wing mounted ft^2
    Embraer.W_fw = 3062; % Weight of feul in the wing lb
    Embraer.q = 0.5*0.00186821*374^2; % dynamic pressure lb/ft^2

%Define Variable Cessna 172
    Cessna172.K_dw = 1; % Wing design factor
    Cessna172.K_vs = 1; % Varible sweep factor
    Cessna172.W_dg = 2200; % Design gross weight (lbs or other units 
    Cessna172.Nz = 5.7; % Ultimate Load factor %ASSUMED
    Cessna172.S_w = 175; % Wing reference area (sq ft or other units)
    Cessna172.A = 7.32; % Aspect ratio
    Cessna172.t_c_root = .12; % Thickness-to-chord ratio at the root
    Cessna172.lambda = 0.75; % Taper ratio
    Cessna172.Lambda = 0; % Sweep angle (in radians)
    %Cessna172.S_csw = 2*11.55; % Control surface area wing mounted ft^2
    Cessna172.W_fw = 252; % Weight of feul in the wing lb
    Cessna172.q = 0.5*0.00207956*249.333^2; % dynamic pressure lb/ft^2

%Define Variable Cessna 180  
    Cessna180.K_dw = 1; % Wing design factor
    Cessna180.K_vs = 1; % Varible sweep factor
    Cessna180.W_dg = 2650; % Design gross weight (lbs or other units 
    Cessna180.Nz = 5.7; % Ultimate Load factor %ASSUMED
    Cessna180.S_w = 175; % Wing reference area (sq ft or other units)
    Cessna180.A = 7.428; % Aspect ratio
    Cessna180.t_c_root = .12; % Thickness-to-chord ratio at the root
    Cessna180.lambda = 0.75; % Taper ratio
    Cessna180.Lambda = 0; % Sweep angle (in radians)
    %Cessna172.S_csw = 2*11.55; % Control surface area wing mounted ft^2
    Cessna180.W_fw = 390; % Weight of feul in the wing lb
    Cessna180.q = 0.5*0.00207956*249.333^2; % dynamic pressure lb/ft^2
%% Wing weight calculation
%W_wing_attack = calculateWingWeightAttack(Embraer.K_dw, Embraer.K_vs, Embraer.W_dg, Embraer.Nz, Embraer.S_w, Embraer.A, Embraer.t_c_root, Embraer.lambda, Embraer.Lambda, Embraer.S_csw);
%
%% Cargo/ Transport
%W_wing_cargo = calculateWingWeightCargo(Embraer.W_dg, Embraer.Nz, Embraer.S_w, Embraer.A, Embraer.t_c_root, Embraer.lambda, Embraer.Lambda, Embraer.S_csw);
% General Aviation
W_wing_GA_EMB = calculateWingWeightGA(Embraer);
W_wing_simple_EMB = 2.5*Embraer.S_w;

W_wing_GA_CSN172 = calculateWingWeightGA(Cessna172);
W_wing_simple_CSN172 = 2.5*Cessna172.S_w;

W_wing_GA_CSN180 = calculateWingWeightGA(Cessna180);
W_wing_simple_CSN180 = 2.5*Cessna180.S_w;
%% Simplified method


%% Fuselage Weight

% Fighter/Attack
%W_fuselage = calculateFuselageWeightAttack(K_dwf, W_dg, Nz, L, D, W);

% Cargo/ Transport


% General Aviation




function W_wing = calculateWingWeightAttack(K_dw, K_vs, W_dg, Nz, S_w, A, t_c_root, lambda, Lambda, S_csw)
    % Function to calculate wing weight based on input parameters
    % Inputs:
    %   K_dw      - Wing design factor
    %   K_vs      - Vertical load factor
    %   W_dg      - Design gross weight
    %   Nz        - Load factor
    %   S_w       - Wing reference area
    %   A         - Aspect ratio
    %   t_c_root  - Thickness-to-chord ratio at the root
    %   lambda    - Taper ratio
    %   Lambda    - Sweep angle at 25% MAC (in radians)
    %   S_csw     - Control surface area (or fraction of wing area)
    %
    % Output:
    %   W_wing - Calculated wing weight

    % Calculate the wing weight using the given formula
    W_wing = 0.0103 * K_dw * K_vs * (W_dg * Nz)^0.5 * S_w^0.622 * A^0.785 * (t_c_root)^-0.4 ...
             * (1 + lambda)^0.05 * (cos(Lambda))^-1.0 * S_csw^0.04;
end
function W_fuselage = calculateFuselageWeightAttack(K_dwf, W_dg, Nz, L, D, W)
    % Function to calculate fuselage weight based on input parameters
    % Inputs:
    %   K_dwf - Fuselage design factor
    %   W_dg  - Design gross weight
    %   Nz    - Ultimate Load factor
    %   L     - Fuselage strctural length
    %   D     - Fuselage strctural depth
    %   W     - Fuselage strctural width
    %
    % Output:
    %   W_fuselage - Calculated fuselage weight

    % Calculate the fuselage weight using the given formula
    W_fuselage = 0.499 * K_dwf * W_dg^0.35 * Nz^0.25 * L^0.5 * D^0.849 * W^0.685;
end
function W_wing = calculateWingWeightCargo(W_dg, Nz, S_w, A, t_c_root, lambda, Lambda, S_csw)
    % Function to calculate wing weight based on input parameters
    % Inputs:
    %   W_dg      - Design gross weight
    %   Nz        - Load factor
    %   S_w       - Wing reference area
    %   A         - Aspect ratio
    %   t_c_root  - Thickness-to-chord ratio at the root
    %   lambda    - Taper ratio
    %   Lambda    - Sweep angle (in radians)
    %   S_csw     - Control surface area (or fraction of wing area)
    %
    % Output:
    %   W_wing - Calculated wing weight

    % Calculate the wing weight using the given formula
    W_wing = 0.0051 * (W_dg * Nz)^0.557 * S_w^0.649 * A^0.5 * (t_c_root)^-0.4 ...
             * (1 + lambda)^0.1 * (cos(Lambda))^-1.0 * S_csw^0.1;
end
function W_wing = calculateWingWeightGA(Plane)
    % Function to calculate wing weight for a general aviation aircraft
    % Inputs:
    %   S_w    - Wing reference area
    %   W_fw   - Weight of the fuel in the wing
    %   A      - Aspect ratio
    %   Lambda - Sweep angle at 25% MAC (in radians)
    %   q      - Dynamic pressure at cruise, lb/ft^2
    %   lambda - Taper ratio
    %   t_c    - Thickness-to-chord ratio
    %   Nz     - Load factor
    %   W_dg   - Design gross weight
    %
    % Output:
    %   W_wing - Calculated wing weight

    % Calculate the wing weight using the given formula
    W_wing = 0.036 * Plane.S_w^0.758 * Plane.W_fw^0.0035 * (Plane.A / cos(Plane.Lambda)^2)^0.6 ...
             * Plane.q^0.006 * Plane.lambda^0.04 * (100 * Plane.t_c_root / cos(Plane.Lambda))^-0.3 ...
             * (Plane.Nz * Plane.W_dg)^0.49;
end