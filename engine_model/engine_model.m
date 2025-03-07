function [percent_power_req] = engine_model(h,V,CL,CD,S)
%main engine model script
%ASSUME MAX RPM IS 2000
%   INPUTS:
        %h: altitude [ft]
        %V: speed in knots
        %S: wing area [ft^2]
%   OUTPUTS:
        %percent power required
%   STEPS:

    %calculated power required
    %check engine + prop model is able to meet that power requirement
    %calculate % power required if yes
    P_max = 1600 * 550; % lbf*ft/s
    rpm_max = 2000;
%calculate power required
    [~, ~, rho, ~] = atmosphere(h);
    V_fts = V * 1.68781; %kts to ft/s
    q = 0.5 * rho * V_fts^2;
    P_req = CD * q * S * V_fts

    %check if engine model can meet power requirement
    torque_percent = interpolateTorque(h, 0);
    [eta, ~] = interpolate_eta(V,rpm_max,h);
    
    T_max = P_max / (rpm_max * (pi/30)); %lbf-ft
    T_avi = (torque_percent/100) * T_max;

    P_avi = T_avi * (rpm_max * (pi/30)) * eta %lbf*ft/s

    if P_avi<P_req
        error('Engine model unable to provide adequate power');
    end

    %return output
    percent_power_req = (P_req/P_avi) * 100;

end

