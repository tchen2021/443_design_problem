function [deltaepsilon_deltaalpha] = downwash(Lambda,AR,lambda, tailfactor, armfactor,h_H,b,l_H)
%DOWNWASH Solves the figures K_A, K_lambda and K_H


A=AR; % Keep notation consistent between equations and our inputs
%% Solve figure K_A
K_A=(1/A)- (1/(1+A^1.7));

%% Solve figure get K_B
K_lambda =(10-(3*lambda))/7;

%% Sovle figure get K_H
K_H= (1-abs(h_H/b))/ (((2*l_H)/b)^(1/3));

%% final calculation
deltaepsilon_deltaalpha = 4.44 * (K_A * K_lambda* K_H * (cosd(Lambda)^0.5))^1.19;

