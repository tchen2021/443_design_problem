%% Mission Profile 1 Weight Estimation
%{

Table of Contents:  

- Section 1: WS and WPL
                takes a range of WS and WPL
- Section 2: Plots
                Carpet plot for WS and WPL
                
- Section 3: Rcr and WPL
                take a range of Rcr and WPL
%}
clear
close all
format short
clc
%% Section 1: WS and WPL

% Inputs
    % Wing Loading
        %WS = linspace(30, 100,36);
        WSLowerLimit = 60;
        WSUpperLimit = 60;
        WSInterval = 1;

            nWS = 1 + (WSUpperLimit-WSLowerLimit)/WSInterval;
            WS = linspace(WSLowerLimit, WSUpperLimit, nWS)';

    % Weapons Payload Weight (Baseline Weight = 500)
        PayloadLowerLimit = 2000;
        PayloadUpperLimit = 3500;
        PayloadInterval = 500;

            nPayload = 1 + (PayloadUpperLimit-PayloadLowerLimit)/PayloadInterval;
            W_PL = linspace(PayloadLowerLimit, PayloadUpperLimit, nPayload);

    % Operational Radius
        Rcr = 300;
    % VFR reserve time (VFRRT)
        VFRRT = 0.5;        % 30 minutes of reserve flight time (at endurance conditions)

%%
% Storage
    WS_WPL.WS = WS;
    WS_WPL.W_PL = W_PL;

    WS_WPL.W_TO = [];
    WS_WPL.W_E = [];
    WS_WPL.W_F = [];

    WS_WPL.EWF = [];
    WS_WPL.Vcr = [];

for i = 1:length(WS)
    for j = 1:length(W_PL)
        k = 4;
        [W_TO, W_E, W_F, EWF, Vcr] = Profile1_propFunction(WS(i), W_PL(j), VFRRT, Rcr, k);

        WS_WPL.W_TO(i,j) = W_TO;
        WS_WPL.W_E(i,j) = W_E;
        WS_WPL.W_F(i,j) = W_F;

        WS_WPL.EWF(i,j) = EWF;
        WS_WPL.Vcr(i,j) = Vcr;
    end
end
%% WS and WPL Plots
%{
% WS and WPL Carpet Plot
[X1,X2] = meshgrid(WS_WPL.W_PL, WS_WPL.WS);
figure;
offset = 0;
nref = 0;

xoff1 =-20;
yoff1 = -200;

xoff2 = 5;
yoff2 = 100;


carpet(X1,X2, WS_WPL.W_TO, 0, nref, 'b', 'r', Linewidth=2)
    carpetlabel(X1,X2, WS_WPL.W_TO, offset, nref, 1, 0, xoff1, yoff1, Color='r', Fontsize=16)
    carpetlabel(X1,X2, WS_WPL.W_TO, offset, nref, 0, 1, xoff2, yoff2, Color='b', Fontsize=16)
    grid minor
    %xlim =
    %ylim = 
    %xlim(xlims)
    %ylim(ylims)
        ax = gca;
        ax.FontSize = 20; 
        ax.YRuler.Exponent =0;

% Labels        
    ylabel("W_T_O [lb]", FontSize=20);

% Constant Values
        strLabel = {'\bfRecon\rm',"SFC = 0.5", "EWF = 0.53", "W/S: 48 lb/sqft", "Loiter Time: 4 hr", "Drag Index: 0.25"};
        text(6200, 19500, strLabel, fontsize=18);

% Lines
        yline(16179, '--', Color = 'k', LineWidth = 1)
        text(8300, 16179-450, 'W_T_O: 16179 lb', FontSize=18)

%% EWF vs L/D

% Inputs
    % Wing Loading
        WSLowerLimit = 40;
        WSUpperLimit = 76;
        WSInterval = 4;

            nWS = 1 + (WSUpperLimit-WSLowerLimit)/WSInterval;
            WS = linspace(WSLowerLimit, WSUpperLimit, nWS)';
    % Drag Index
        k = [5 4 3 2 1];    % best drag values to worst drag values
        extradrag = 0.0000;
        LDcr_4script = linspace(9.17199982/(1+extradrag), 17.3167353/(1+extradrag), 5);            
    % WPL
        W_PL = 2500;
    % Operational Radius
        Rcr = 300;
    % VFR reserve time (VFRRT)
        VFRRT = 0.5;        % 30 minutes of reserve flight time (at endurance conditions)
    % 
% Storage
    EWF_LD.WS = WS;
    EWF_LD.W_PL = W_PL;

    EWF_LD.W_TO = [];
    EWF_LD.W_E = [];
    EWF_LD.W_F = [];

    EWF_LD.EWF = [];
    EWF_LD.Vcr = [];

for i = 1:length(WS)
    for j = 1:length(k)
        [W_TO, W_E, W_F, EWF, Vcr] = Profile1_propFunction(WS(i), W_PL, VFRRT, Rcr, k(j));

        EWF_LD.W_TO(i,j) = W_TO;
        EWF_LD.W_E(i,j) = W_E;
        EWF_LD.W_f(i,j) = W_F;

        EWF_LD.EWF(i,j) = EWF;
        EWF_LD.Vcr(i,j) = Vcr;
        EWF_LD.LDcr(i,j) = LDcr_4script(k(j));
    end
end

%EWF and LD Carpet Plot
X1 = EWF_LD.EWF;
X2 = EWF_LD.LDcr;
figure; grid on;
offset = 0; nref = 0;
carpet(X1,X2, EWF_LD.W_TO, 0, nref, 'b', 'r')
    carpetlabel(X1',X2', EWF_LD.W_TO', offset, nref, 1, 0, 0, 0)
    carpetlabel(X1',X2', EWF_LD.W_TO', offset, nref, 0, 1, 0, 0)

%}
%% WPL and Rcr
% 
% % Inputs
%     % Wing Loading
%         WS = 66;
%     % Weapons Payload Weight (Baseline Weight = 500)
%         PayloadLowerLimit = 2500;
%         PayloadUpperLimit = 3500;
%         PayloadInterval = 100;
% 
%             nPayload = 1 + (PayloadUpperLimit-PayloadLowerLimit)/PayloadInterval;
%             W_PL = linspace(PayloadLowerLimit, PayloadUpperLimit, nPayload);
%     % Operational Radius
%         Rcr = 300;
%     % VFR reserve time (VFRRT)
%         VFRRT = 0.5;        % 30 minutes of reserve flight time (at endurance conditions)
% 
% % Storage
%     Rcr_WPL.WS = WS;
%     Rcr_WPL.W_PL = W_PL;
% 
%     Rcr_WPL.W_TO = [];
%     Rcr_WPL.W_E = [];
%     Rcr_WPL.W_F = [];
% 
%     Rcr_WPL.EWF = [];
%     Rcr_WPL.Vcr = [];
% 
% for i = 1:length(Rcr)
%     for j = 1:length(W_PL)
%         [W_TO, W_E, W_F, EWF, Vcr] = Profile1_propFunction(WS, W_PL(j), VFRRT, Rcr(i));
% 
%         Rcr_WPL.W_TO(i,j) = W_TO;
%         Rcr_WPL.W_E(i,j) = W_E;
%         Rcr_WPL.W_f(i,j) = W_F;
% 
%         Rcr_WPL.EWF(i,j) = EWF;
%         Rcr_WPL.Vcr(i,j) = Vcr;
%     end
% end

%% 
% 
% 
% %%Old Plots
% 
% figure
% hold on; grid minor
% plot(WS,1600./W_TO, 'b', LineWidth=2)
% plot(WS,1600./W_TO, 'o', Color='b')
% title("W_t_o vs W/S")
% xlabel("W/S [lb/ft^2]")
% ylabel("T/W [lb]")
% figure;
% hold on; grid minor
% plot(WS,EWF, 'b', LineWidth=2)
% plot(WS,EWF, 'o', Color='b')
% title("EWF vs W/S")
% xlabel("W/S [lb/ft^2]")
% ylabel("EWF [lb]")
% WTO_vs_WS_turboprop_profile1 = [WS' W_TO'];
% 
% 
% 
% % Power Loading Plot
    % P = 1600;
    % figure; hold on; grid on;
    % plot(WS,P./WS_WPL.W_TO(:,1))
    % plot(WS,P./WS_WPL.W_TO(:,2))
    % plot(WS,P./WS_WPL.W_TO(:,3))
    % plot(WS,P./WS_WPL.W_TO(:,4))
    % legend("2000lb", "2500lb", "3000lb", "3500lb")
    % title("Power/Weight vs WS")

