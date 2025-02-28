


A=[1,2,5,10,40]; % L/D

B=[.6,.7,.8,.9,1]; %Drag (C_D)

answer=polyfit(A,B,6);

polyval(answer,.6)