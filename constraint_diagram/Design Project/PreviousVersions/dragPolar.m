function [CD,Cl,Cl_maxL,Cl_max] = dragPolar(W,rho,S,V,V_max,V_s,AR,e,CD_0)
%Calculates the drag polar of aircraft given inputs

% Calculate max values given max speeds
Cl_max= (2*W)/(rho*S*V_max^2);
Cl_maxL= (2*W)/(rho*S*V_s^2); 

Cl=(2*W)./(rho*S.*V.^2);
CD= CD_0 + Cl.^2./(pi*AR*e);
end