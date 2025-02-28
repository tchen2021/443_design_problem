function SFC = turbopropSFC(SFC0, V, h, dT)
    % turbopropSFC calculates the SFC of a turboprop engine with temperature offset
    % Inputs:
    %   SFC0 - Reference SFC at sea level static conditions (lbm/hp/hr)
    %   V    - Airspeed in knots
    %   h    - Altitude in feet
    %   dT   - Temperature offset from ISA (in °F)
    %
    % Output:
    %   SFC  - Specific Fuel Consumption (lbm/hp/hr)

    % Constants
    V_ref = 280; % Reference airspeed in knots
    rho0 = 0.002377; % Sea-level air density [slugs/ft^3]
    C1 = 0.05; % Airspeed effect coefficient
    C2 = 0.2; % Altitude effect coefficient

    % Calculate temperature at altitude (ISA) in Rankine with offset
    T0 = 518.67; % Sea level temperature in Rankine
    L = 0.00356616; % Temperature lapse rate in Rankine/ft
    T = T0 - L * h + (dT * 1.8); % Apply temperature offset (°F to °R)

    % Calculate pressure at altitude (ISA) in psf
    P0 = 2116.22; % Sea-level pressure [psf]
    P = P0 * (1 - (L * h) / T0)^5.2561; % Pressure at altitude [psf]

    % Calculate air density (slug/ft^3)
    R = 1576; % Specific gas constant for air [ft*lbf/(slug*R)]
    rho = P / (R * T);
    sqrt_theta = T0/T;
    gamma = 1.4;
    a = sqrt(gamma*R*T);

    M = V/a;

    % Calculate SFC
%     SFC = SFC0 * (1 + C1 * (V / V_ref) + C2 * (rho0 / rho));
    SFC = (0.3+0.9*M)*sqrt_theta;
end
