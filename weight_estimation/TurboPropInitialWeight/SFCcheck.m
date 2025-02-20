
clear
clc
close all
format short


SFC0 = 0.5;
V = [260 280 300];
h = linspace(0,25000, 1001);
dT = 0;
for j = 1:length(V)
    V_input = V(j);
    for i = 1:length(h)
        h_input = h(i);
        SFC = turbopropSFC(SFC0, V_input, h_input, dT)
        SFC_plot(i,j) = SFC;
    end
end
figure; grid minor; hold on;
plot(SFC_plot(:,1), h, LineWidth=2)
plot(SFC_plot(:,2), h, LineWidth=2)
plot(SFC_plot(:,3), h, LineWidth=2)

ylabel("Altitude, ft")
xlabel("SFC")
ax = gca;
ax.YRuler.Exponent = 0;
legend("260 kts", "280 kts", "300 kts")








    % turbopropSFC calculates the SFC of a turboprop engine with temperature offset
    % Inputs:
    %   SFC0 - Reference SFC at sea level static conditions (lbm/hp/hr)
    %   V    - Airspeed in knots
    %   h    - Altitude in feet
    %   dT   - Temperature offset from ISA (in °F)