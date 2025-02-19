% 2 profile overlayed plot

clear
close all
format short
clc


data1 = load("Profile1_2500_40_76.mat")
data2 = load("Profile2_3500_40_76.mat")

EWF_LD_1 = data1.EWF_LD;
EWF_LD_2 = data2.EWF_LD;

X1_1 = EWF_LD_1.EWF;
X2_1 = EWF_LD_2.LDcr;

%[X1,X2] = meshgrid(EWF_LD_1.EWF, EWF_LD_1.LDcr);
figure; hold on;
offset = 0.0;
nref = 0;

xoff1 =0.0;
yoff1 = -3500;

xoff2 = 0.065;
yoff2 = 7000;


carpet(-0.000+X1_1,X2_1, EWF_LD_1.W_TO, 0.00015, nref, 'r', 'b', Linewidth=1.5, Linestyle = '-')

carpet(X1_1,X2_1, EWF_LD_2.W_TO, 0.00015, nref, 'r', 'b', Linewidth=1.5, Linestyle = '--')

    %carpetlabel(X1_1',X2_1', EWF_LD_2.W_TO', offset, nref, 1, 0, xoff1, yoff1, Color='r', Fontsize=16)
    %carpetlabel(X1_1',X2_1', EWF_LD_1.W_TO', offset, nref, 0, 1, xoff2, yoff2, Color='b', Fontsize=16)
    grid minor
    xlims = [0.496 0.58];   %[0.52 0.54];
    ylims = [10000 34000];
        yticks(10000:2000:34000);
    xlim(xlims)
    ylim(ylims)
        ax = gca;
        ax.FontSize = 20; 
        ax.YRuler.Exponent =0;
set(groot, 'DefaultAxesFontName', 'Calibri');   % Change axes font
set(groot, 'DefaultTextFontName', 'Calibri');   % Change text font
% Labels        
    ylabel("W_T_O [lb]", FontSize=20);

% % Constant Values
%         strLabel = {'\bfRecon\rm',"SFC = 0.5", "EWF = 0.53", "W/S: 48 lb/sqft", "Loiter Time: 4 hr", "Drag Index: 0.25"};
%         text(6200, 19500, strLabel, fontsize=18);
% 
% % Lines
%         yline(16179, '--', Color = 'k', LineWidth = 1)
%         text(8300, 16179-450, 'W_T_O: 16179 lb', FontSize=18)