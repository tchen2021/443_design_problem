function [TW_TO] = takeoffDistance(sigma,Cl_maxTO,S_TOFL,TOP25)
%This function calculates takeoff distance as a ratio of Thrust/Weight
    % This function gets a value of T/W so that it can be plotted in a
    % constraint diagram

%% Calculate Take off weight
TW_TO= @(WS) (WS/ (sigma.*Cl_maxTO*TOP25*S_TOFL) );
end