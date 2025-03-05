 

%%proffesors altitude function
function [p, T, rho, a] = atmospheresendi(h)
% atmosphere(h) function compute p, T, rho and a, giving the altitude h
% p : pressure in lb/ft^2
% T : temperature in degree Rankine
% rho : density in slug/ft^3
% a : speed of sound
% h : altitude in ft, h<=104990ft
%% Constants
Cv = 4290; % constant volume specific heat for the air
R = 1716; % ideal gas constant
Cp = R+Cv; % constant pressure specific heat for air
gamma = Cp/Cv; % specific heat ratio
T_sl = 518.69; % sea level temperature
p_sl = 2116.2; % sea level pressure
%% temperature ratio
if h<=36089
theta = 1-(6.875*10^-6)*h;
elseif h>36089 && h<=65617
theta = 0.75189;
elseif h>65617 && h<=104990
theta = 0.75189+(1.0577*(10^-6))*(h-65617);
else
theta=nan;
warning('altitude exceeded')
end
%% pressure density
delta = ( 1-(6.875*(10^-6)).*h ).^5.2561.* ((0 <= h) & (h <= 36089)) +...
(0.2234*exp(4.806*(10^-5).*(36089-h))).* ((36089 < h )& (h <= 65617)) +...
(3.174716*(10^-6).*(0.75189+1.0577.*(10^-6).*(h-65617)).^-34.164).*((65617<h) & (h<=104990));
%% p, T and rho
T = T_sl.*theta;
p = p_sl.*delta;
rho = p./(R.*T);
a = sqrt(gamma*R.*T);
end