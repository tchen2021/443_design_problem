function [eta, D] = interpolate_eta(V,n,h)
%INTERPOLATES 4-BLADED PROPELLER DATA TO FIND eta GIVEN J = V/nD
%   INPUTS
    %V = freestream speed [kts]
    %n = max rpm of engine (2000rpm) [ft]
    %h = altitude [ft]
%   OUTPUTS
    %eta: propeller efficiency
    %D:   propeller diameter. Fixed value once solved for high altitude
    %cruising speed of 280kts at 25000ft

    
    
    propdata = importfileGeneral("prop_data.csv", 2); %import digitized data
    J_bounds = [propdata(1,1), propdata(end,1)];      %find the domain of data
    
    %calculate J (advance ratio)
    V = V * 1.68781; %kts to ft/s
    [~, ~, ~, a] = atmosphere(h);
    proptip_speedlimit = 0.8 * a; %to avoid superonic effects on prop
    % rpm_max = (proptip_speedlimit/(D/2))/(pi/30);
    %D = (2*proptip_speedlimit)/(n * (pi/30)); %use rad/s here for rotational speed
    D = 7.8;        %ft
    J = V/((n*(1/60)*D));                     %use rotations/s for rotational speed

    %check if advance ratio is out of bounds
    if J< J_bounds(1) || J > J_bounds(2)
        error('Advance Ratio out of bounds');
    end
    
    p = polyfit(propdata(:,1),propdata(:,2),15);
    
    %J = 0.2;
    eta = polyval(p, J);
end

