clear
clc
close all
set(groot, 'DefaultAxesFontName', 'Calibri');   % Change axes font
set(groot, 'DefaultTextFontName', 'Calibri');   % Change text font

dt = 15;
disp("Turboprop")
% Example 1: Calculate SFC at sea level and 200 knots
sfc0 = 0.5;
sfc1 = turbopropSFC(sfc0, 280, 300,dt);
fprintf('SFC at sea level and 200 knots: %.4f lbm/hp/hr\n', sfc1);

% Example 2: Calculate SFC at 20,000 ft and 300 knots
sfc2 = turbopropSFC(sfc0, 300, 20000,dt);
fprintf('SFC at 20,000 ft and 300 knots: %.4f lbm/hp/hr\n', sfc2);

% Example 3: Calculate SFC at 10,000 ft and 250 knots
sfc3 = turbopropSFC(sfc0, 250, 10000,dt);
fprintf('SFC at 10,000 ft and 250 knots: %.4f lbm/hp/hr\n', sfc3);

disp("turbojet")
tfsc0 = 0.85;
% Example 1: Calculate TSFC at sea level and 300 knots
tsfc1 = turbojetTSFC(tfsc0, 300, 0,dt);
fprintf('TSFC at sea level and 300 knots: %.4f lbm/lbf/hr\n', tsfc1);

% Example 2: Calculate TSFC at 25,000 ft and 500 knots
tsfc2 = turbojetTSFC(tfsc0, 500, 25000,dt);
fprintf('TSFC at 25,000 ft and 500 knots: %.4f lbm/lbf/hr\n', tsfc2);

% Example 3: Calculate TSFC at 15,000 ft and 400 knots
tsfc3 = turbojetTSFC(tfsc0, 400, 15000,dt);
fprintf('TSFC at 15,000 ft and 400 knots: %.4f lbm/lbf/hr\n', tsfc3);

n = 6;
V_prop = 280; %[kts]
V_jet = 310;  %[kts]
V_array = linspace(0,350,n);
h = linspace(0,30000,n); %[ft]
SFC_array = zeros(n,n);
TSFC_array = zeros(n,1);
for i =1:n
    for j = 1:n
    SFC_array(i,j) = turbopropSFC(sfc0,V_array(j),h(i),dt);
    end
end

[X1, X2] = meshgrid(h,V_array);
fig = figure('Visible', 'off');  % No figure window appears
figure;
hold on
grid on
offset = -3;
nref = 0;

carpet(X1,X2, SFC_array, offset, nref, 'r', 'b', Linewidth = 2);
%carpetlabel(X1,X2, SFC_array, offset, nref, 'r', 'b');
% ylim([0.55 0.85])
xlim([-5000 34000])
ax = gca;
ax.YRuler.Exponent = 0;
ax.FontSize = 20;
hold off
width = 16;  % width in inches
height = 9;  % height in inches
set(gcf, 'Units', 'inches', 'Position', [1 1 width height]);

% Adjust paper size for saving
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height], 'PaperSize', [width height]);

% Save the figure as an SVG file
saveas(gcf, 'sfc_plot.svg');  % Alternatively: exportgraphics(gcf, 'example_plot.svg')

% figure
% hold on
% title("Prop")
% ylabel("Alt.")
% xlabel("SFC")
% plot(SFC_array,h,LineWidth=1)
% 
% grid on
% hold off

% figure
% hold on
% title("Jet")
% ylabel("Alt.")
% xlabel("SFC")
% plot(TSFC_array,h,LineWidth=1)
% grid on
% hold off
