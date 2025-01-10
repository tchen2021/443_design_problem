function[MTOWL,ClMaxList] = landingDistance(S_l)
%Outputs a vector of maximum takeoff wing loadings 
% with corresponding ClMaxes

%Takes inputs of:
%Minimum Field Length
%Operating at ISA+10ºC @MSL
T0 = 288.15;%K
%deltaT = 10;
delT = 10;
rhoSTD = 2.376892 * 10^-3;
rho = rhoSTD/(1+(delT/T0));

% ClMaxList = linspace(1.9,3.3,3);
ClMaxList = 2.6; % CHOSEN CL
wingLoadingRatio = .8; %Landing at 80% MTOW
MTOWL = ClMaxList;
for i = 1:length(ClMaxList)
    MTOWL(i) = ((S_l*rho*ClMaxList(i))/(2*.3))/wingLoadingRatio;
end