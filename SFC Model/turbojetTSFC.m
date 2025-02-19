function TSFC = turbojetTSFC(TSFC0, V, h, dT)
    % turbojetTSFC calculates the TSFC of a turbojet engine with temperature offset
    % Inputs:
    %   TSFC0 - Reference TSFC at sea level static conditions (lbm/lbf/hr)
    %   V     - Airspeed in knots
    %   h     - Altitude in feet
    %   dT    - Temperature offset from ISA (in °F)
    %
    % Output:
    %   TSFC  - Thrust Specific Fuel Consumption (lbm/lbf/hr)

    % Constants
    V_ref = 400; % Reference airspeed in knots
    rho0 = 0.002377; % Sea-level air density [slugs/ft^3]
    C1 = 0.03; % Airspeed effect coefficient
    C2 = 0.15; % Altitude effect coefficient

    % Calculate temperature at altitude (ISA) in Rankine with offset
    T0 = 518.67; % Sea level temperature in Rankine
    L = 0.00356616; % Temperature lapse rate in Rankine/ft
    T = T0 - L * h + (dT * 1.8); % Apply temperature offset (°F to °R)

    % Calculate pressure at altitude (ISA) in psf
    P0 = 2116.22; % Sea-level pressure [psf]
    P = P0 * (1 - (L * h) / T0)^5.2561; % Pressure at altitude [psf]

    % Calculate air density (slug/ft^3)
    R = 1716; % Specific gas constant for air [ft*lbf/(slug*R)]
    rho = P / (R * T);

    % Calculate TSFC
    TSFC = TSFC0 * (1 + C1 * (V / V_ref) + C2 * (rho0 / rho));
end