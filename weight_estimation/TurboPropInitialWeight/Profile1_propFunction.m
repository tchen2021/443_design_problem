%% Profile1_propFunction
%{
Inputs: WS      Wing Loadings
        W_PL    Weapons Payload Weight
        VFRRT   VFR Reserve Time
        Rcr     Operational Radius

Outputs: W_TO   Takeoff Weight
         W_E    Empty Weight
         W_F    Fuel Weight
         EWF    Empty Weight Fraction
         Vcr    Cruise Velocity

Units:  Weight      lb
        Distance    nautical miles
        Speed       kts
        Time        hours
        SFC         lb / hp*hr
%}
function [W_TO, W_E, W_F, EWF, Vcr] = Profile1_propFunction(WS, W_PL, VFRRT, Rcr,k)
%% Aerodynamic Values
        extradrag = 0.0000;
    % Cruise L/D
        LDcr = linspace(10.2/(1+extradrag), 16.1/(1+extradrag), 5);     % bad (1) to good (0)

    % Loiter L/D
        LDlt = linspace( 8.8/(1+extradrag), 14/(1+extradrag), 5);

    % Loiter Velocity
        rho_lt = 0.001066;       % slugs/ft^3, 25000ft   
        rho_cr = 0.001756; %0.001066;       % slugs/ft^3, 25000ft
        
        Cl_lt_0 = 1;     % CL_Endurance_0;                     % Cl loiter at lowest drag configuration
        Cl_lt_1 = 1;     % CL_Endurance_1;                     % Cl loiter at highest drag configuration
        Vlt_0 = sqrt( (2*WS) / (rho_lt*Cl_lt_0) );  % ft/s
        Vlt_1 = sqrt( (2*WS) / (rho_lt*Cl_lt_1) );  % ft/s

        Vlt = (linspace(Vlt_0,Vlt_1, length(LDlt)))* 0.592484  % knots

        Vcr = 0.592484 * sqrt((2*WS) / (rho_cr*0.4));       % knots0.378788

%% Mission Profile 1: Full Operational Radius + Full Loiter Time (Constant Payload) 
 
% Weight Ratios                
    W1TO = 0.990;
    W21 = 0.990;
    W32 = 0.990;
    W43 = 0.980;
%     W54 = 1 / (exp((R_cr * SFC_cr)/(LD_cr * etap * 325)));
%     W65 = 1 / (exp((E * V_lt * SFC_lt)/(LD_lt * etap * 325)));
%     W76 = 1 / (exp((R_cr * SFC_cr)/(LD_cr * etap * 325)));
    W87 = 0.990;
    W98 = 1;        % alternate replaced extra loiter time; 1 / (exp((R_alt * SFC_alt)/(LD_alt * etap * 325)));
    W109 = 0.995;
 
    W_ratios = W1TO*W21*W32*W43*W87*W98*W109;   % combine constant weight ratios for simplicity
    
    % Operational Radius [nm]
        Rcr = Rcr;
    % Loiter Time [hrs]
        E = 4;
%% Engine & Propeller Variables
    % Cruise & Loiter Specific Fuel Consumption
        SFCcr = 0.5; SFClt = SFCcr;
    % Propeller Efficiency
        etap = 0.85;
%% EWF as a function of WS
   EWF =  -0.030 + 1.009*WS^-0.145;

%% System of Equations
    % Drag Polar
        DragPolarResolution = length(LDcr); % 
% Structure for each calculated weight
    % Rows    -> Variable 1
    % Cols    -> Variable 2
    % "Pages" -> Sensitivy Analysis Variable
    % (rows, cols, pages)
        Q2_WS_check.W_to = [];
        Q2_WS_check.W_e = [];
        Q2_WS_check.W_f = [];

%%for DragPolarIndex = 1:DragPolarResolution
    % Fsolve options
        options = optimoptions('fsolve','Display','off');
    % Initial Guess
        initial_guess = [15000, 7800, 3380];
    % Set Drag Index
        %k = 4;

    for j = 1:length(Rcr )         % Calculate weights at every range
        for  i= 1:length(W_PL)                      % Calculate weights at every W_PL

            MFF = W_ratios * ((1 / (exp((Rcr(j) * SFCcr)          /(LDcr(k) * etap * 325)))) * ...
                              (1 / (exp(((E+VFRRT) * Vlt(k) * SFClt)/(LDlt(k) * etap * 325)))) * ...
                              (1 / (exp((Rcr(j) * SFCcr)          /(LDcr(k) * etap * 325)))));

            equations = @(x) [ -x(1) + x(2) + x(3) + W_PL(i) + 500; % 500 for pilots and sensors
                               -x(1) + (1/EWF)*x(2);
                               x(1)*(1 - MFF) - x(3)];
    
            solutions = fsolve(equations, initial_guess, options);

            Q2_WS_check.W_to(j, i) = solutions(1);  % Takeoff Weight
            Q2_WS_check.W_e(j, i) = solutions(2);  % Empty Weight
            Q2_WS_check.W_f(j, i) = solutions(3);   % Fuel Weight
        end
    end

% Results
    W_TO = Q2_WS_check.W_to(j, i);
    W_E = Q2_WS_check.W_e(j, i);
    W_F = Q2_WS_check.W_f(j, i);

%end

%DragIndexPlot = flip([0 0.25 0.5 0.75 1]);
%[X1, X2] = meshgrid (W_pl, R_cr);




%% Plots
% figure; hold on; grid on;
% ylims = [4000 22000];
% xlims = [6000 17000];
% 
% carpet(X1, X2, WPL_Drag.W_to(:,:), 35, 0, 'r', 'b', Linewidth = 2);
%     xlim(xlims)
%     ylim(ylims)
%         ax = gca;
%         ax.FontSize = 20; 
%         ax.YRuler.Exponent =0;
%         grid on;
% 
% % Labels        
%     ylabel("W_T_O [lb]", FontSize=20);
% 
%     text(15500, 12200, "W_P_L [lb]", color = 'b', FontSize=20);
%     carpetlabel(X1, X2, WPL_Drag.W_to(:,:), 0, 0, 1, 0, 11400, 0, 'Color', 'b', 'FontSize', 18)
% 
%     text(12800, 20400, "Range [nm]", color = 'r', FontSize=20)
%     carpetlabel(X1,X2, WPL_Drag.W_to(:,:), 35, 1, 0, 1, -100, 400, 'Color', 'r', 'FontSize', 18)
% 
% % Constant Values
%         strLabel = {'\bfRecon\rm',"SFC = 0.5", "EWF = 0.53", "W/S: 48 lb/sqft", "Loiter Time: 4 hr", "Drag Index: 0.25"};
%         text(6200, 19500, strLabel, fontsize=18);
% 
% % Lines
%         yline(16179, '--', Color = 'k', LineWidth = 1)
%         text(8300, 16179-450, 'W_T_O: 16179 lb', FontSize=18)

