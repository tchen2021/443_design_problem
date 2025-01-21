%% Class 2 estimation

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
% Define the variables
K_dw = 1; % Wing design factor
K_vs = 1; % Varible sweep factor
W_dg = 12500; % Design gross weight (lbs or other units as needed)
Nz = 3.75; % Ultimate Load factor
S_w = 313; % Wing reference area (sq ft or other units)
A = 7.98; % Aspect ratio
t_c_root = .12; % Thickness-to-chord ratio at the root
lambda = 1; % Taper ratio
Lambda = 0; % Sweep angle (in radians)
S_csw = 200; % Control surface area wing mounted ft^2

% Wing weight calculation
W_wing_attack = calculateWingWeightAttack(K_dw, K_vs, W_dg, Nz, S_w, A, t_c_root, lambda, Lambda, S_csw);

% Cargo/ Transport
W_wing_cargo = calculateWingWeightCargo(W_dg, Nz, S_w, A, t_c_root, lambda, Lambda, S_csw);
% General Aviation
W_wing_GA = calculateWingWeightGA(S_w, W_fw, A, Lambda, q, lambda, t_c, Nz, W_dg);

%% Fuselage Weight

% Fighter/Attack
W_fuselage = calculateFuselageWeightAttack(K_dwf, W_dg, Nz, L, D, W);

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
function W_wing = calculateWingWeightGA(S_w, W_fw, A, Lambda, q, lambda, t_c, Nz, W_dg)
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
    W_wing = 0.036 * S_w^0.758 * W_fw^0.0035 * (A / cos(Lambda)^2)^0.6 ...
             * q^0.006 * lambda^0.04 * (100 * t_c / cos(Lambda))^-0.3 ...
             * (Nz * W_dg)^0.49;
end