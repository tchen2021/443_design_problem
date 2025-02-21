function [M] = machCalc(V,h)
%Solves for the mach number at the given speed

% std atmo
[p, T, rho, a] = atmosphere(h);
% Constants for air
gamma=1.4;
R =1717;
V=knots2ftperS(V);
% Solve for Mach number
M= V/sqrt(T*R*gamma); 
end % end of function