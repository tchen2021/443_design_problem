%% Longitudinal Aerodynamic Deviatives Calculator 
%Authors: Tyson Chen & Roman Niemiec
%This main script takes geometric data of the plane design pulled from an excel sheet named
% "DragBuildUp.xls and calculates aerodynamic peformance of it using the drag buildup
% notes on Canvas
clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                          DATA SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adding derivatives functions and data
%addpath derivatives;
addpath(genpath('drag'))
addpath lift
addpath downwash
addpath zeroliftalpha
addpath lift_slope
addpath XFoil_Results


% Import data file
%DragBuildUp = importfileTEST('DragBuildUp.xlsx');
DragBuildUp = importfileTEST('DragBuildUp.xlsx');
% Take individual variables from excel sheet
for i=1:height(DragBuildUp)
    assignin('base',DragBuildUp{i,2},DragBuildUp{i,3});
end
                    
j=1;
%%% ITERATOR for CL CD 
airfoil = 'NACA 2412';
airfoil_tailhorizontal = 'NACA 0012';       %importer tool does not import strings
airfoil_tailvertical = 'NACA 0012';

%manually adjusting the alpha iteration ranges
begin = -10;
ending = 20;

%initialize data structures to store drag values 
C_L_struct.winglift = zeros(1, (ending-begin)+1);



for alpha=begin:1:ending 

% Get variables from Speed and altitude
[p, T, rho, a] = atmosphere(h); % get temp in [R]
[M] = machCalc(V,h); % Mach number based on speed of air outside airplane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          LIFT SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift slope - Wing
% Solve Graph B11A for K value
[K_Wing,Re] = B11A(h,V,C_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Wing] = B11B(t_root,t_tail,C_root,C_tip);

%%% 5. Calculate beta
beta_Wing=sqrt(1-M^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Wing=(1.05/beta_Wing)*K_Wing*Cl_alpha_theory_Wing;

%%% 6. Calculate k
k_Wing=(beta_Wing*Cl_alpha_Wing)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);% xaxis=(A/k_Wing)*sqrt(beta_Wing^2+tand(Lambda_half)^2)
% [CL_alpha_Wing] = B12(xaxis)/A;
% CL_alpha_Wing= ( (2*pi)/(2+sqrt( ((A^2 * beta_Wing^2)/2)*(1+ (tand(Lambda_half)^2)/beta_Wing^2 )+4))) *A; % Function mentioned on graph 
CL_alpha_Wing=  ((A*2*pi) / (2+ sqrt( ((A^2)*(beta_Wing^2))/k_Wing^2 *(1+(tand(Lambda_half)^2)/beta_Wing^2)+4)));



%% Lift slope - TAIL SECTION %%
% Solve Graph B11A for K value
[K_Horiz,Re] = B11A(h,V,C_Vert_root);

% Solve Graph B11B for Cl_alpha_theory
% Calc wing thickness ratio
[Cl_alpha_theory_Horiz] = B11B(t_Vert_root,t_Vert_tip,C_Vert_root,C_Vert_tip);

%%% 5. Calculate beta
beta_Horiz=sqrt(1-M^2);

%%% 4. Calculate Cl_alpha
Cl_alpha_Horiz=(1.05/beta_Horiz)*K_Horiz*Cl_alpha_theory_Horiz;

%%% 6. Calculate k
k_Horiz=(beta_Horiz*Cl_alpha_Horiz)/(2*pi);

%%% 7. Find CL_alpha
% xaxis=(A/K_Wing)*sqrt(beta_Wing^2 + tand(Lambda_half)^2);
% [CL_alpha_Wing] = B12(xaxis)/A;
CL_alpha_Horiz= ( (2*pi)/(2+sqrt( ((A^2 * beta_Horiz^2)/2)*(1+ (tand(Lambda_Tail_half)^2)/beta_Horiz^2 )+4))) *A_T; % Function mentioned on graph 

%% zero-lift angle of attack
% alpha_0W_root = deg2rad(-5); %zero lift angle of attack, 2D, at the wing root
% iprime_r = deg2rad(-3);      %geometric incidence of the wing at the root
%J = ;             %empirical factor

alpha_0W = zeroliftalpha(iprime_r, epsilon, alpha_0W_root, AR, lambda);

%% Downwash
delepsilondelalpha = downwash(Lambda_quarter, AR, lambda, tailfactor, armfactor,h_H,b,l_H); % how the wake effects of the wing affect the flow of the tail

%% Fuselage Lift
%largest fuse area 9786.86
% 2786.86




% Solve Lift value for lift body
C_LB = fuselift(lf_D, alpha, alpha_0B, S, D, S_P_x0, M, la, lf);

%%% Wing body lift %%%
a_W = CL_alpha_Wing;
C_LWB = wingbodylift(C_LB, a_W, alpha, alpha_0W, Se_S, D_b); %alpha in rad

%%% horizontal tail lift %%%
a1 = CL_alpha_Horiz;          %lift slope of the tail
C_LT = a1 .* ((deg2rad(alpha)*(1-delepsilondelalpha)) - (epsilon_0 + i_T))* (S_T/S); % Calculate with alpha in rads

%total lift
C_L_struct.winglift(j) = C_LWB-C_LB;
C_L(j) = C_LWB + C_LT;
% end % end of test for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          DRAG SECTION                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wing Drag
% For loop to iterate through 4 sections of our wing
S_i=[S1,S2,S3,S4];
cbarbari=[cbarbar1,cbarbar2,cbarbar3,cbarbar4];

C_D0W(j) = iterator(airfoil,alpha,iprime_r, epsilon,cbarbari, S_i,V,rho,Se,S,h);
[C_DiW(j), C_DiW_delta1(j), C_DiW_delta2(j)] = C_DiWfun(lambda, Lambda_quarter, C_LWB-C_LB, AR);

%Tail
alpha_T = alpha + i_T - delepsilondelalpha .* alpha;   %the angle of attack of horizontal tail
cbarbari_horizontal = [cbarbar1_horizontal cbarbar2_horizontal cbarbar3_horizontal cbarbar4_horizontal];
S_i_horizontal = [S1_horizontal S2_horizontal S3_horizontal S4_horizontal];

cbarbari_vertical = [cbarbar1_vertical cbarbar2_vertical cbarbar3_vertical cbarbar4_vertical];
S_i_vertical = [S1_vertical S2_vertical S3_vertical S4_vertical];
% C_D0T = 2.* ((sum(C_DTi .* S_Ti)) / (sum(S_Ti))) * (S_Te/S);   %subscript T refers to the tail (horizontal or vertical)
C_D0T(j)= iterator(airfoil_tailhorizontal, alpha, i_T, epsilon_0, cbarbari_horizontal, S_i_horizontal, V,rho, STe, S, h); % drag caused by horizontal tail WE AINT DOING THIS YET
C_D0TV(j)= iterator(airfoil_tailvertical, 0, 0, 0, cbarbari_vertical, S_i_vertical, V,rho, STe_V, S, h); % drag caused by vertical tail WE AINT DOING THIS YET
%deltaC_D0T =               %additional drag due to deflection of the
%elevator not accounted for now


% delta3_f1 = ;          %interpolated from figure flap 1
% delta3_f2 = ;          %interpolated from figure flat 2
% 
% delta_C_D0T = delta1 * delta2 * (delta3_f1 - delta3_2); %additional drag caused by a plain flap (elevator)


C_DiT(j) = C_DiTfun(C_LT, A_T, S, S_T, t_c, beta_Horiz, cf_C); % UNCOMMENT SOON


%% Fuselage and Nacelles
C_D0F(j) = C_D0Ffun(Re, S, S_A_wet, S_B_wet, S_C_wet, S_CAB, 't', lf, D, Lambda_f);       %assumes smooth surface turbulent flow
C_DiB(j) = C_LB*(deg2rad(alpha)-(deg2rad(alpha_0B)));          %induced drag coefficient of the fuselage w/ alpha in radians

C_D(j)=C_DiB(j)+C_DiT(j)+C_D0T(j)+C_D0TV(j)+C_DiW(j)+C_D0W(j) + C_D0F(j);

%% Interference Drag

%wing-fuselage
deltaC_DWf(j) = 0.05 * (C_D0F(j) + C_DiB(j)); 

%wing-nacelle
%n = ;   %number of nacelles
%deltaC_Dwn = 0.05 * (C_D0N + C_DiN) * n;
deltaC_Dwn(j) = 20e-4;             %clean configurations 20 drag counts
%Tail-Fuselage
% n_c = ;     %number of corners between tail and fuselage
% c_j = ;     %chord of the tail at the junction
% deltaC_Dtf = n_c * (0.8 * (t/c)^3 - 0.0005) * (c_j^2/S);
deltaC_Dtf(j) = 10e-4;                %streamlined tail design 10 drag counts

%deltaC_Dtt = (n_c/2) * (17*(t/c)^4 - 0.05*(t/c)^2) * (c_j^2/S); %for two
%tails

%cooling and air intake
% mdot = ;        %mass flow rate of air inside the intake
% q = ;           %dynamic pressure
% deltaV =  ;     %variation of the velocity between the intake and the outlet
% C_D_cooling = (mdot*deltaV)/(q*S);
C_D_cooling(j) = 30e-4;      %air intake and cooling 30 drag counts


% Update the CD with interferance drag
C_D_Interferance(j) = deltaC_DWf(j)+deltaC_Dwn(j)+deltaC_Dtf(j)+C_D_cooling(j);




%%% TOTAL DRAG 
C_D(j)= C_D(j)+C_D_Interferance(j);
j=j+1;
end % end of for loop through

[~,CD0_index] =min(abs(C_L));  %returns the index for the CD0 value for later use

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%                                                          Landing gear, external armaments, flaps                                                                       %%%%%%%%%%
%%%%%%%%%%                                                                                                                                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Munitions
%%%% ASSUMING GBU12 500 LB BOMB (actually weights 510 lbs)
C_D_GBU12=C_D_bombWithPylonNOINTDRAG *(area_GBU12/S); % Calc drag with respect to our wing size to proportion the drag counts
C_D_WeaponsAttack= C_D_GBU12*bombNumAtt ;
C_D_WeaponsRecon=C_D_GBU12*bombNumRec;

C_D_250Tank=C_D_dropTank*(area_250DropTank/S);
C_D_75Tank=C_D_dropTankWing*(area_75DropTank/S);
C_D_120Tank=C_D_dropTank*(area120DropTank/S);
C_D_175Tank=C_D_dropTank*(area175DropTank/S);

C_D_FuelAttack= 0;
% C_D_FuelRecon= 0;
C_D_FuelRecon= C_D_175Tank;

% ISR Camera
C_D_Camera=C_D_ISR*(area_ISR/S);

% Total extra drag from external stuff - TYSON FOR THE ATTACK AND RECON
% MISSION JUST ADD THESE ON!!!!
C_D_Attack=C_D_WeaponsAttack+C_D_FuelAttack+C_D_Camera;
C_D_Recon=C_D_WeaponsRecon+C_D_FuelRecon+C_D_Camera;


%% calculating drag contributions as % of total drag
totDrag_clean = C_D(CD0_index);
totDrag_attack = totDrag_clean + C_D_Attack;
totDrag_recon = totDrag_clean + C_D_Recon;

fprintf('Total drag counts: \nClean: %d\nAttack: %d\nRecon: %d', round(totDrag_clean*10^4), round(totDrag_attack*10^4), round(totDrag_recon*10^4));

%Wing
dragContribution.clean.Wing = (C_D0W(CD0_index) + C_DiW(CD0_index)) / totDrag_clean;
dragContribution.attack.Wing = (C_D0W(CD0_index) + C_DiW(CD0_index)) / totDrag_attack;
dragContribution.recon.Wing = (C_D0W(CD0_index) + C_DiW(CD0_index)) / totDrag_recon;

%Fuselage
dragContribution.clean.Fuselage = (C_DiB(CD0_index) + C_D0F(CD0_index)) / totDrag_clean;
dragContribution.attack.Fuselage = (C_DiB(CD0_index) + C_D0F(CD0_index)) / totDrag_attack;
dragContribution.recon.Fuselage = (C_DiB(CD0_index) + C_D0F(CD0_index)) / totDrag_recon;

%Tail
dragContribution.clean.Tail = (C_DiT(CD0_index)+C_D0T(CD0_index)+C_D0TV(CD0_index)) / totDrag_clean;
dragContribution.attack.Tail = (C_DiT(CD0_index)+C_D0T(CD0_index)+C_D0TV(CD0_index)) / totDrag_attack;
dragContribution.recon.Tail = (C_DiT(CD0_index)+C_D0T(CD0_index)+C_D0TV(CD0_index)) / totDrag_recon;

%Interference
dragContribution.clean.Interference = C_D_Interferance(CD0_index) / totDrag_clean;
dragContribution.attack.Interference = C_D_Interferance(CD0_index) / totDrag_attack;
dragContribution.recon.Interference = C_D_Interferance(CD0_index) / totDrag_recon;

%External Payloads
dragContribution.clean.Payload = 0 / totDrag_clean;
dragContribution.attack.Payload = C_D_Attack / totDrag_attack;
dragContribution.recon.Payload = C_D_Recon / totDrag_recon;

%Displaying Results
fprintf('\n\nTotal drag counts (CLEAN): \nWing: %0.4f\nFuselage: %0.4f\nTail: %0.4f\nInterference: %0.4f\nExternal Payloads: %0.4f\n', dragContribution.clean.Wing, (dragContribution.clean.Fuselage), (dragContribution.clean.Tail), (dragContribution.clean.Interference), dragContribution.clean.Payload)
fprintf('\n\nTotal drag counts (ATTACK): \nWing: %0.4f\nFuselage: %0.4f\nTail: %0.4f\nInterference: %0.4f\nExternal Payloads: %0.4f\n', dragContribution.attack.Wing, (dragContribution.attack.Fuselage), (dragContribution.attack.Tail), (dragContribution.attack.Interference), dragContribution.attack.Payload)
fprintf('\n\nTotal drag counts (RECON): \nWing: %0.4f\nFuselage: %0.4f\nTail: %0.4f\nInterference: %0.4f\nExternal Payloads: %0.4f\n', dragContribution.recon.Wing, (dragContribution.recon.Fuselage), (dragContribution.recon.Tail), (dragContribution.recon.Interference), dragContribution.recon.Payload)


%% saving particular results into .mat file for convenience
%save("MainFunctionData.mat", "C_L_struct");

%% finding max range and endurance CL, CD
%drag buildup: C_L, C_D    1x13
%competitive analysis: CL, CD 100x4

load old_drag_polar.mat

%find the max L/D CL point for drag buildup
[placeholder, idx_max(1)] = max(C_L(:)./C_D(:));
CL_max(1) = C_L(idx_max(1));
CD_max(1) = C_D(idx_max(1));
LD_max(1) = CL_max(1)/CD_max(1);

%find the max L/D CL point for competitive analysis clean configuration
[placeholder, idx_max(2)] = max(CL(:,1)./CDtot(:,1));
CL_max(2) = CL(idx_max(2),1);
CD_max(2) = CDtot(idx_max(2),1);
LD_max(2) = CL_max(2)/CD_max(2);

%find the max L/D CL point for competitive analysis for the high drag case
[placeholder, idx_max(3)] = max(CL(:,2)./CDtot(:,2));
CL_max(3) = CL(idx_max(3),2);
CD_max(3) = CDtot(idx_max(3),2);
LD_max(3) = CL_max(3)/CD_max(3);

%find the max L/D CL point for competitive analysis for the 0.25 drag index
[placeholder, idx_max(4)] = max(CL_third(:)./CDtot_third(:));
CL_max(4) = CL_third(idx_max(4));
CD_max(4) = CDtot_third(idx_max(4));
LD_max(4) = CL_max(4)/CD_max(4);

%find the max L/D CL point for attack mission drag buildup (limiting case
%between recon and attack)
[placeholder, idx_max(5)] = max(C_L(:)./ (C_D(:)+C_D_Attack));
CL_max(5) = C_L(idx_max(5));
CD_max(5) = C_D(idx_max(5)) + C_D_Attack;
LD_max(5) = CL_max(5)/CD_max(5);


%find the max endurance for prop drag buildup
[placeholder, idx_max_endurance(1)] = max((C_L(:)).^1.5 ./ C_D(:));
CL_endurance(1) = C_L(idx_max_endurance(1));
CD_endurance(1) = C_D(idx_max_endurance(1));
LD_endurance(1) = CL_endurance(1)/CD_endurance(1);


%find the max endurance for prop competitive analysis clean configuration
[placeholder, idx_max_endurance(2)] = max((CL(:,1)).^1.5 ./ CDtot(:,1));
CL_endurance(2) = CL(idx_max_endurance(2),1);
CD_endurance(2) = CDtot(idx_max_endurance(2),1);
LD_endurance(2) = CL_endurance(2)/CD_endurance(2);

%find the max endurance for prop competitive analysis for the high drag
%case
[placeholder, idx_max_endurance(3)] = max((CL(:,2)).^1.5 ./ CDtot(:,2));
CL_endurance(3) = CL(idx_max_endurance(3),2);
CD_endurance(3) = CDtot(idx_max_endurance(3),2);
LD_endurance(3) = CL_endurance(3)/CD_endurance(3);

%find the endurance L/D CL point for competitive analysis for the 0.25 drag index
[placeholder, idx_max_endurance(4)] = max(CL_third(:).^1.5./CDtot_third(:));
CL_endurance(4) = CL_third(idx_max_endurance(4));
CD_endurance(4) = CDtot_third(idx_max_endurance(4));
LD_endurance(4) = CL_endurance(4)/CD_endurance(4);

%find the endurance L/D CL point for attack mission drag buildup (limiting case
%between recon and attack)
[placeholder, idx_max_endurance(5)] = max(C_L(:).^1.5./(C_D(:) + C_D_Attack));
CL_endurance(5) = C_L(idx_max_endurance(5));
CD_endurance(5) = C_D(idx_max_endurance(5)) + C_D_Attack;
LD_endurance(5) = CL_endurance(5)/CD_endurance(5);




%% plotting

blue = '#2E5F7F';
red = '#A03232';
green = '#33692F';
orange = '#D97843';
aqua = '#499CD0';
cp_gold = '#BD8B13';
cp_green = '#154734';
purple = '#ff00ff';

% Check the data
% Plot CL vs Alpha
% plot(begin:ending,C_L); xlabel("Alpha");ylabel("CL");grid on;

polarWidth = 4;
sz = 500;
linethickness = 4;
color = 'magenta';
color2 = [0.4660 0.6740 0.1880];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT THE MUTHAFUCKIN DRAG POLA
plot(C_D,C_L, '--', 'Color',green, 'LineWidth', 4); hold on; 
set(groot, 'DefaultAxesFontName', 'Calibri');   % Change axes font
set(groot, 'DefaultTextFontName', 'Calibri');   % Change text font
%scatter(CD_max(1), CL_max(1), sz, 'square',  'MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max range point
%scatter(CD_endurance(1),CL_endurance(1), sz, 'o','MarkerEdgeColor', purple,  'LineWidth', linethickness) %plot max endurance point
%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlabel("C_D");ylabel("C_L");grid on;
%plot formatting
fontSize_axes = 26;
fontSize_text = 28;
fontSize_subtitles = 28;
offset = 0.03; %offset from the horizontal line
yLineWidth = 3;

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Good case (competitve analysis)
plot(CDtot(:,1), CL(:,1), 'Color', red, 'LineWidth', 6) %left polar
%scatter(CD_max(2), CL_max(2), sz, 'square',  'MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max range point
%scatter(CD_endurance(2),CL_endurance(2), sz, 'o','MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max endurance point
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bad case (competitive analysis)
%hold on;
plot(CDtot(:,2), CL(:,2), 'Color', red, 'LineWidth', 6) %right polar
%scatter(CD_max(3), CL_max(3), sz,  'square',  'MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max range point
%scatter(CD_endurance(3),CL_endurance(3), sz, 'o','MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max endurance point
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIRTY CONFIG DRAG POLARS (ATTACK AND RECON)
plot(C_D+C_D_Attack, C_L, ':', 'Color', orange, 'LineWidth', polarWidth) % Attack Mission
plot(C_D+C_D_Recon, C_L, '--', 'Color', aqua, 'LineWidth', polarWidth) % Recon Mission
%scatter(CD_max(5), CL_max(5), sz,  'square',  'MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max range point
%scatter(CD_endurance(5),CL_endurance(5), sz, 'o','MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max endurance point
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%drag index 0.25 polar

%plot(CDtot_third, CL_third,'--', 'Color', orange, 'LineWidth', polarWidth) 
%scatter(CD_max(4), CL_max(4), sz,  'square',  'MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max range point
%scatter(CD_endurance(4),CL_endurance(4), sz, 'o','MarkerEdgeColor', purple, 'LineWidth', linethickness) %plot max endurance point
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%labeling the polars
%text2 = sprintf('Highest Drag Index = 1\n{C_D}_0 = %.3f\ne = %.2f\nA = %.2f', CD0(2), e(2), A(2)); %for highet drag case

%text(0.1, 0.525, text2, 'FontName', 'Abolition', 'FontSize', fontSize_subtitles)

%text1 = sprintf('Clean Drag Index = 0\n{C_D}_0 = %.3f\ne = %.2f\nA = %.2f', CD0(1), e(1), A(1)); %for clean conditions

%text(0.0025, 0.67, text1, 'FontName', 'Abolition', 'FontSize', fontSize_subtitles)

%plot formatting

xlim([0 0.1]);           % X-axis limit from 0 to 0.2
ylim([0 1.2]);           % Y-axis limit from 0 to 1.4

% Set major ticks for grid lines
xticks(0:0.02:0.2);      % Major ticks every 0.02 on X-axis
yticks(0:0.2:1.4);       % Major ticks every 0.2 on Y-axis

% Enable grid only at major ticks
grid on;                 % Enable grid lines at major ticks
ax = gca;                % Get current axes
ax.GridLineStyle = '-';  % Solid line for major grid

% Set minor ticks without grid lines
ax.XMinorTick = 'on';           % Enable minor ticks on X-axis
ax.YMinorTick = 'on';           % Enable minor ticks on Y-axis
ax.MinorGridLineStyle = 'none'; % Turn off minor grid lines

% Define minor tick intervals
ax.XAxis.MinorTickValues = 0:0.005:0.2;  % Minor ticks every 0.005 on X-axis
ax.YAxis.MinorTickValues = 0:0.05:1.4;   % Minor ticks every 0.05 on Y-axis

% Label axes
xlabel('C_D');           % Label for X-axis
ylabel('C_L', Rotation= 0);           % Label for Y-axis

% Set font size and other formatting adjustments to match style
ax.FontSize = fontSize_axes;        % Adjust font size
ax.LineWidth = 1;        % Set axis line width

%legend('Current','Preliminary (Optimistic)');


%% plotting the lift curve slope
figure

scatter(deg2rad(begin:1:ending), C_L_struct.winglift);
