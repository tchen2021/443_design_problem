clear
close all
format short
clc

WS = linspace(1,100);
W_TO = [];
EWF = [];
for i = 1:length(WS)
    [W_TO(i), EWF(i)] = Q2_WS_Check(WS(i));
end  
%%
figure
hold on; grid minor
plot(WS,W_TO, 'b', LineWidth=2)
plot(WS,W_TO, 'o', Color='b')
title("W_t_o vs W/S")
xlabel("W/S [lb/ft^2]")
ylabel("W_T_O [lb]")




%%
figure
hold on; grid minor
plot(WS,1600./W_TO, 'b', LineWidth=2)
plot(WS,1600./W_TO, 'o', Color='b')
title("W_t_o vs W/S")
xlabel("W/S [lb/ft^2]")
ylabel("T/ [lb]")
figure;
hold on; grid minor
plot(WS,EWF, 'b', LineWidth=2)
plot(WS,EWF, 'o', Color='b')
title("EWF vs W/S")
xlabel("W/S [lb/ft^2]")
ylabel("EWF [lb]")
WTO_vs_WS_turboprop_profile1 = [WS' W_TO'];
