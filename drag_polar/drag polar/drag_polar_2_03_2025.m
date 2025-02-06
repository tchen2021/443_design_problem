%% DRAG POLAR 
%housekeeping
close all 
clear 
clc

%% load data and L/D from competitive asssessment
dragpolar_data = load('dragpolar_data.mat');
CL_market = dragpolar_data.ComparisionTable2S3.CL; %CL data
CD_market = dragpolar_data.ComparisionTable2S3.CD; %CD data
Market_names = dragpolar_data.ComparisionTable2S3.Name;
% i think the range should be about 8 to 13 - from sam 



%% first drag polar (props)
% Example point 
%CL_dat = 0.40;
%CD_dat = 0.025 + 0.01;


% LD(1) = 12.15224; % T-6B 
LD(1) = 16; % 13.3675 A-37 Dragonfly jet max L/D to verify jets

x(:,1) = linspace(0, 0.1);
y(:,1) = LD(1).*x(:,1); % this is just y = mx + b; to plot the line 


% recall the equation for a drag polar 
% CDtot = CD0 + 1/k * CL^2 

% pick your range of CL, at this point its just a guess 
CL(:,1) = linspace(0, 1.5); 
CD0(1) = 0.018; % also just a guess 0.0182

e(1) = 0.85; % guess with this one too to get it close  %0.3
A(1) = 7; %T-6B 6.515662651
k(1) = (pi * e(1) * A(1)); 

CDtot(:,1) = CD0(1) + 1/k(1) * CL(:,1).^2; 


%find the max L/D CL point
[placeholder, idx_max(1)] = max(CL(:,1)./CDtot(:,1))
CL_max(1) = CL(idx_max(1),1);
CD_max(1) = CDtot(idx_max(1),1);
LD_max(1) = CL_max(1)/CD_max(1);

%find the carson speed for prop
[placeholder, idx_carson_speed(1)] = max((CL(:,1)).^0.5 ./ CDtot(:,1))
CL_carson(1) = CL(idx_carson_speed(1),1);
CD_carson(1) = CDtot(idx_carson_speed(1),1);
LD_carson(1) = CL_carson(1)/CD_carson(1);

%find the max range for prop
[placeholder, idx_range(1)] = max((CL(:,1)) ./ CDtot(:,1))
CL_range(1) = CL(idx_range(1),1);
CD_range(1) = CDtot(idx_range(1),1);
LD_range(1) = CL_range(1)/CD_range(1);

%find the max endurance for prop
[placeholder, idx_max_endurance(1)] = max((CL(:,1)).^1.5 ./ CDtot(:,1))
CL_endurance(1) = CL(idx_max_endurance(1),1);
CD_endurance(1) = CDtot(idx_max_endurance(1),1);
LD_endurance(1) = CL_endurance(1)/CD_endurance(1);
%% second drag polar (props)

% lets make the bad line 
% e_bad = 0.27;  % iscold plays with this value so that the lift values are the same across the horiz lines
% A_bad = 5.2; 
% kbad = (pi * e_bad * A_bad);
% CD0_bad = 0.035; 
% CDtot_bad = CD0_bad + 1/kbad * CL(:,1).^2; 
% [LD_bad, idx] = max(CL(:,1)./CDtot_bad);

% LD(2) = 9.722; % OV-10D 
LD(2) = 9.7; % %OV-10 Bronco 11.6 L-39 Albotross

x(:,2) = linspace(0, 0.15);
y(:,2) = LD(2).*x(:,2); % this is just y = mx + b; to plot the line 


% recall the equation for a drag polar 
% CDtot = CD0 + 1/k * CL^2 

% pick your range of CL, at this point its just a guess 
CL(:,2) = linspace(0, 1.5); 
CD0(2) = 0.028; % also just a guess 0.028

e(2) = 0.7; % guess with this one too to get it close  0.21
A(2) = 5.35; %OV-10 (defunct) 
k(2) = (pi * e(2) * A(2)); 


CDtot(:,2) = CD0(2) + 1/k(2) * CL(:,2).^2; 

%find the max L/D CL point
[placeholder, idx_max(2)] = max(CL(:,2)./CDtot(:,2))
CL_max(2) = CL(idx_max(2),2);
CD_max(2) = CDtot(idx_max(2),2);
LD_max(2) = CL_max(2)/CD_max(2);

%find the carson speed for prop
[placeholder, idx_carson_speed(2)] = max((CL(:,2)).^0.5 ./ CDtot(:,2))
CL_carson(2) = CL(idx_carson_speed(2),2);
CD_carson(2) = CDtot(idx_carson_speed(2),2);
LD_carson(2) = CL_carson(2)/CD_carson(2);

%find the max range for prop
[placeholder, idx_range(2)] = max((CL(:,2)) ./ CDtot(:,2))
CL_range(2) = CL(idx_range(2),2);
CD_range(2) = CDtot(idx_range(2),2);
LD_range(2) = CL_range(2)/CD_range(2);

%find the max endurance for prop
[placeholder, idx_max_endurance(2)] = max((CL(:,2)).^1.5 ./ CDtot(:,2))
CL_endurance(2) = CL(idx_max_endurance(2),2);
CD_endurance(2) = CDtot(idx_max_endurance(2),2);
LD_endurance(2) = CL_endurance(2)/CD_endurance(2); 


%% third drag polar index 0.25


% recall the equation for a drag polar 
% CDtot = CD0 + 1/k * CL^2 

% pick your range of CL, at this point its just a guess 
CL_third = linspace(0, 1.5)'; 
% CD0_third = CDtot(1) + 0.25*(CDtot(2)-CDtot(1)); % also just a guess 0.028

% e_third = 0.5; % guess with this one too to get it close  0.21
% A_third = 6.5; %OV-10 (defunct) 
% k_third = (pi * e_third * A_third); 


%CDtot_third = CD0_third + 1/k_third * CL_third.^2; 
CDtot_third = CDtot(:,1) + 0.25.*(CDtot(:,2)-CDtot(:,1)); % also just a guess 0.028

%%
%inputting color hex codes
blue = '#2E5F7F';
red = '#A03232';
green = '#33692F';
orange = '#D97843';
aqua = '#499CD0';
cp_gold = '#BD8B13';
cp_green = '#154734';
purple = '#ff00ff';


%% plotting props (unlabeled planes, no jets)
figure 
hold on

%plot formatting
fontSize_axes = 26;
fontSize_text = 28;
fontSize_subtitles = 28;
offset = 0.03; %offset from the horizontal line
yLineWidth = 3;
polarWidth = 4;

%plot(x(:,1), y(:,1), 'b', 'LineWidth', 2) %max L/D line for left polar
%text(0.051, 0.65, 'T-6B Best L/D ','FontSize', fontSize_subtitles, 'FontName', 'Calibri', 'Rotation', 55, Color='k')
plot(CDtot(:,1), CL(:,1), 'Color', red, 'LineWidth', polarWidth) %left polar
plot(CDtot_third, CL_third, '--', 'Color', 'k', 'LineWidth', polarWidth) %0.25 drag index
scatter(CDtot(idx_max(1),1), CL(idx_max(1),1), 120,  "o",'LineWidth', 1.5, MarkerEdgeColor = orange)

%plot sample cruise condition points from competitive assessment
for i = 1:11
    if i== 9 || i==10 || i == 8 || i == 3 %skip AT-802A and PC-9 due to insufficient data and AU-24 (outlier)
        continue

    elseif i == 1 || i == 2 %increase marker size and unique color for main competitor AT-6 and A29
        scatter(CD_market(i), CL_market(i), 120, 'filled', MarkerFaceColor = purple, MarkerEdgeColor = purple) %purple #ff00ff
        %text(CD_market(i) + 0.005, CL_market(i), Market_names(i), 'FontName', 'Calibri', 'FontSize', 12, 'Color', '#4B4B4B'); % Adjust 0.
    else

        scatter(CD_market(i), CL_market(i), 90, 'filled', MarkerFaceColor = '#4B4B4B', MarkerEdgeColor = '#4B4B4B') %dark gray
        %text(CD_market(i) + 0.005, CL_market(i), Market_names(i), 'FontName', 'Calibri', 'FontSize', 12, 'Color', '#4B4B4B'); % Adjust 0.02 as needed

    end
end

% %plotting jets scat
% for i = 12:18
%     scatter(CD_market(i), CL_market(i), 'filled', MarkerFaceColor = '#00FF00', MarkerEdgeColor = '#00FF00') %dark gray
%     %text(CD_market(i) + 0.005, CL_market(i), Market_names(i), 'FontName', 'Calibri', 'FontSize', 12, 'Color', '#4B4B4B'); % Adjust 0.02 as needed
% 
% end

%plot(x(:,2), y(:,2), "Color", '#77AC30', 'LineWidth', 2) %max L/D line for right polar
%text(0.073, 0.65, 'OV-10D Best L/D ','FontSize', fontSize_subtitles, 'Rotation', 48, Color='k')    
plot(CDtot(:,2), CL(:,2), 'Color', blue, 'LineWidth', polarWidth) %right polar
scatter(CDtot(idx_max(2),2), CL(idx_max(2),2), 120, "o",'LineWidth', 1.5, MarkerEdgeColor = orange) %"#D95319"

%plot best range line
yline(CL(idx_max(1),1), '--', 'Color', green, 'LineWidth', yLineWidth)
%text(0.10, CL(idx_max(1),1)+offset, 'Best Range', 'FontName', 'Calibri','FontSize', fontSize_text, Color='k')
%text(0.08, CL(idx_max(1),1)+offset, sprintf('L/D = %.1f', LD_max(1)), 'FontName', 'Calibri', 'FontSize', fontSize_text, Color='k')  
%text(0.08, CL(idx_max(1),1)-offset, sprintf('L/D = %.1f', LD_max(2)), 'FontName', 'Calibri','FontSize', fontSize_text,  Color='k') 

%plot max endurance line 
yline(CL(idx_max_endurance(1),1), '--', 'Color', green, 'LineWidth', yLineWidth)
%text(0.06, CL(idx_max_endurance(1),1)+offset, 'Max Endurance', 'FontName', 'Calibri','FontSize', fontSize_text, Color='k')
scatter(CDtot(idx_max_endurance(1),1), CL(idx_max_endurance(1),1), 120, 'o', 'LineWidth', 1.5, MarkerEdgeColor = orange)
scatter(CDtot(idx_max_endurance(2),2), CL(idx_max_endurance(2),2), 120, 'o','LineWidth', 1.5, MarkerEdgeColor = orange)
%text(0.04, CL(idx_max_endurance(1),1)+offset, sprintf('L/D = %.1f', LD_endurance(1)), 'FontName', 'Calibri', 'FontSize', fontSize_text, Color='k')  
%text(0.04, CL(idx_max_endurance(1),1)-offset, sprintf('L/D = %.1f', LD_endurance(2)), 'FontName', 'Calibri', 'FontSize', fontSize_text, Color='k') 



%labeling the polars
text2 = sprintf('Highest Drag Index = 1\n{C_D}_0 = %.3f\ne = %.2f\nA = %.2f', CD0(2), e(2), A(2)); %for highet drag case

%text(0.1, 0.525, text2, 'FontName', 'Calibri', 'FontSize', fontSize_subtitles)

text1 = sprintf('Clean Drag Index = 0\n{C_D}_0 = %.3f\ne = %.2f\nA = %.2f', CD0(1), e(1), A(1)); %for clean conditions

%text(0.0025, 0.67, text1, 'FontName', 'Calibri', 'FontSize', fontSize_subtitles)

% %plotting max speed CL
% WS = 73; %lb/ft^2
% q = 0.5*0.002378*(280*1.68781)^2; 
% CL_max_speed(1) = WS/q;
% yline(CL_max_speed(1), '--', 'Color', green, 'LineWidth', yLineWidth)
% text(0.17, CL_max_speed(1)+offset, 'Max Speed', 'FontName', 'Calibri','FontSize', fontSize_text, Color='k')
% idx_max_speed(1) = 19;
% scatter(CDtot(idx_max_speed(1),1), CL(idx_max_speed(1),1), 120, 'o', 'LineWidth', 1.5, MarkerEdgeColor = orange);
% scatter(CDtot(idx_max_speed(1),2), CL(idx_max_speed(1),2), 120, 'o', 'LineWidth', 1.5, MarkerEdgeColor = orange);
% LD_max_speed(1) = CL(idx_max_speed(1),1) / CDtot(idx_max_speed(1),1); %calculate LD at min. cruise condition
% LD_max_speed(2) = CL(idx_max_speed(1),2) / CDtot(idx_max_speed(1),2);
% text(0.08, CL(idx_max_speed(1),1)+offset, sprintf('L/D = %.1f', LD_max_speed(1)), 'FontName', 'Calibri', 'FontSize', fontSize_text, Color='k')
% text(0.08, CL(idx_max_speed(1),2)-offset, sprintf('L/D = %.1f', LD_max_speed(2)), 'FontName', 'Calibri', 'FontSize', fontSize_text, Color='k')

%plot formatting

xlim([0 0.13]);           % X-axis limit from 0 to 0.2
ylim([0 1.1]);           % Y-axis limit from 0 to 1.4

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


set(gca, 'FontName', 'Calibri')