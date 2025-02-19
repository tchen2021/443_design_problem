% Example 1: Calculate SFC at sea level and 200 knots
sfc1 = turbopropSFC(0.5, 200, 0);
fprintf('SFC at sea level and 200 knots: %.4f lbm/hp/hr\n', sfc1);

% Example 2: Calculate SFC at 20,000 ft and 300 knots
sfc2 = turbopropSFC(0.5, 300, 20000);
fprintf('SFC at 20,000 ft and 300 knots: %.4f lbm/hp/hr\n', sfc2);

% Example 3: Calculate SFC at 10,000 ft and 250 knots
sfc3 = turbopropSFC(0.5, 250, 10000);
fprintf('SFC at 10,000 ft and 250 knots: %.4f lbm/hp/hr\n', sfc3);

function SFC = turbopropSFC(SFC0, V, h)
    % turbopropSFC calculates the SFC of a turboprop engine in lbm/hp/hr
    % Inputs:
    %   SFC0 - Reference SFC at sea level static conditions (lbm/hp/hr)
    %   V    - Airspeed in knots
    %   h    - Altitude in feet
    %
    % Output:
    %   SFC  - Specific Fuel Consumption (lbm/hp/hr)
    
    % Constants
    V_ref = 300; % Reference airspeed in knots
    rho0 = 0.002377; % Sea-level air density [slugs/ft^3]
    C1 = 0.05; % Airspeed effect coefficient
    C2 = 0.2; % Altitude effect coefficient

    % Calculate temperature at altitude (ISA) in Rankine
    T0 = 518.67; % Sea level temperature in Rankine
    L = 0.00356616; % Temperature lapse rate in Rankine/ft
    T = T0 - L * h; % Temperature in Rankine

    % Calculate pressure at altitude (ISA) in psf
    P0 = 2116.22; % Sea-level pressure [psf]
    P = P0 * (1 - (L*h)/T0)^5.2561; % Pressure at altitude [psf]

    % Calculate air density (slug/ft^3)
    R = 1716; % Specific gas constant for air [ft*lbf/(slug*R)]
    rho = P / (R * T);

    % Calculate SFC
    SFC = SFC0 * (1 + C1 * (V / V_ref) + C2 * (rho0 / rho));
end
