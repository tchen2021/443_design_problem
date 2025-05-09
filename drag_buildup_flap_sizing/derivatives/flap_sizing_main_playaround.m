%% Flap Sizing Program
%% Tyson Chen
%%%%%%%%%%%%%%%%%%%%%%%%%%
%performs flap sizing based on Young, A.D., "The Aerodynamic Characteristics of Flaps", Aeronautical Research Council Reports and
%Memoranda, R&M 2622, 1953.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%first input a target CL, CL generated by wing, flap span ratio b_f/b
%program varies flap deflection to find chord ratio 
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASSUMPTIONS
%critical case for flap sizing is landing approach, where CL_max is
%needed

%THIS CAN BE REFINED WITH A TAKEOFF MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%UPDATES TO BE MADE
%automatically obtain CL generated by wing based on speed, altitude, alpha
%by using MainFunction.m and developing takeoff model
%%%%%%%%%%%%%%%%%%%%%%%%%%

%housekeeping 
clc; clear; close all;

%adding interpolated .csv from filepath
%addpath flap_digitization;
addpath(genpath('digitizer_files'))

%importing constants from excel sheet
% Import data file
DragBuildUp = importfileTEST('DragBuildUp.xlsx');

%importing CL data to get deltaCL
load drag_polar.mat;

% Take individual variables from excel sheet
for i=1:11      %Plane property section
    assignin('base',DragBuildUp{i,2},DragBuildUp{i,3});
end
%A = 7;
assignin('base',DragBuildUp{15,2},DragBuildUp{15,3}); %import A, aspect ratio
t_c = 0.15;             %tail thickness ratio
lambda = 0.472972973;   %taper ratio
%assignin('base',DragBuildUp{86,2},DragBuildUp{866,3}); %import A, aspect ratio

%inputting required values

%stall speed 76.43 kts
%approach speed 148.42 ft/s 
% V_A = 148.41;                    %in ft/s
% q = 0.5 * 0.002378 * V_A^2;    %approach dynamic pressure at sea level
% CL_W = W/(q*S);                 %input CL generated by the wing at conditions
CL_target = 2.0;                %input targeted CL for wing and flaps combination
CL_max = 1.6;                   %theory of wing sections, p.478 appendix for NACA 2412

%deltaCL = CL_target - 0.9*C_L_struct.winglift(end);
deltaCL = CL_target - 0.9 * CL_max;
%deltaCL = .65;
bf_b = 0.75;                    %input ratio of flap span to wing span
%bf_b = 1;


%% analysis section (interpolating data, solving for chord ratio
flap_type = 'p'; %plain flap
%flap_type = 'f'; %fowler flap
%flap_type = "sl"; %slotted flap
%flap_type = "sp"; %split flap

deflection_start = 0;       %degrees
deflection_end = 50;        %degrees
deflection_interval = 5;   %degrees

deflections = deflection_start:deflection_interval:deflection_end;
init_arr = zeros(length(deflections),1);

[lambda1, lambda2, lambda22, lambda3, delta1, delta2, delta3, K, mu1, mu2, cf_c, deltaCD0, deltaCDi, deltaCM] = deal(init_arr);
FA_F6 = FA_F6fun(A);

                      %iterator index
for j=1:length(deflections)

    lambda2(j) = lambda2fun(deflections(j), t_c, flap_type);
    lambda3(j) = lambda3fun(bf_b);
    %delta1(j) = delta1fun
    delta2(j) = delta2fun(deflections(j), t_c, flap_type);
    %delta3(j) = delta3fun(bf_b);
    %mu1(j) = mu1fun
    mu2(j) = mu2fun(bf_b);
    K(j) = Kfun(bf_b);
    
    %calculate cf_c flap chord ratio
    lambda1(j) = deltaCL / (FA_F6 * lambda2(j) * lambda3(j));
    cf_c(j) = 0.1753;
    delta1(j) = delta1fun(cf_c(j), t_c, flap_type);
    deltaCD0(j) = (delta1(j) * delta2(j)) * (delta3fun(bf_b, lambda)-delta3fun(0, lambda));
    deltaCDi(j) = K(j) * deltaCL^2 / (pi*A);
    deltaCM(j) = -mu1fun(cf_c(j), t_c) * (mu2fun(j)-mu2fun(0)) * deltaCL; 
end

%% Plotting the drag polars accounting for flaps
%importing drag polar from drag buildup DI=0.25 (CDtot_third and CL_third)
load dragpolar_025_edwards.mat

%define CL ranges for takeoff and landing
CL_landing_range = [1.26 2.0]; %estimates for STOL LTA
CL_takeoff_range = [0.26 1.8];
CD0_landing_gear = 0.0075;    %Roman's calculation

%plotting takeoff drag polar
target_index = 9; %flap deflection of 40deg, cf/c=0.17
deflection_landing = deflections(target_index);
deflection_takeoff = 19;
deltaCD0_landing = deltaCD0(target_index);
deltaCDi_landing = deltaCDi(target_index);
deltaCD_landing = deltaCD0_landing + deltaCDi_landing;

deltaCD0_takeoff = deltaCD0_landing * deflection_takeoff/deflection_landing;
deltaCDi_takeoff = deltaCDi_landing * (0.9)^2;
deltaCD_takeoff = deltaCDi_takeoff + deltaCD0_takeoff;

%CD_landing = CDtot_third + deltaCD_landing + CD0_landing_gear;
%CL_landing = CL_third;

%CD_takeoff = CDtot_third + deltaCD_takeoff + CD0_landing_gear;
%CL_takeoff = CL_third; 

% %finding elements in the takeoff range
% idx_takeoff = CL_takeoff >= CL_takeoff_range(1) & CL_takeoff <= CL_takeoff_range(2);  % Logical index of values within range
% CL_takeoff_subset = CL_takeoff(idx_takeoff);          % x values in the range
% CD_takeoff_subset = CD_takeoff(idx_takeoff);          % Corresponding y values
% 
% %finding elements in the landing range
% idx_landing = CL_landing >= CL_landing_range(1) & CL_landing <= CL_landing_range(2);  % Logical index of values within range
% CL_landing_subset = CL_landing(idx_landing);          % x values in the range
% CD_landing_subset = CD_landing(idx_landing);          % Corresponding y values

%extrapolating the CL and CD polars to account for flaps
drag_polar_clean = polyfit(C_L, C_D, 2); %drag polar of clean

CL_takeoff = linspace(1.26, CL_takeoff_range(2));
CD_takeoff = polyval(drag_polar_clean,CL_takeoff) + deltaCD_takeoff + CD0_landing_gear;
K_takeoff_arr = polyfit(CL_takeoff,CD_takeoff,2);
K_takeoff = K_takeoff_arr(1);

CL_landing = linspace(1.26, CL_landing_range(2));
CD_landing = polyval(drag_polar_clean, CL_landing) + deltaCD_landing + CD0_landing_gear;
K_landing_arr = polyfit(CL_landing,CD_landing,2);
K_landing = K_landing_arr(1);

%plotting
figure; hold on;

blue = '#2E5F7F';
red = '#A03232';
green = '#33692F';
orange = '#D87941 ';
aqua = '#499CD0';
cp_gold = '#BD8B13';
cp_olive = '#B38E33';
cp_green = '154734';
purple = '#ff00ff';

polarWidth = 4;
sz = 500;
linethickness = 4;

grid on;
%plot formatting
fontSize_axes = 26;
fontSize_text = 28;
fontSize_subtitles = 28;
offset = 0.03; %offset from the horizontal line
yLineWidth = 3;

plot(C_D, C_L, 'LineWidth', polarWidth);
plot(CD_takeoff, CL_takeoff, 'LineWidth', polarWidth);
plot(CD_landing, CL_landing, 'LineWidth', polarWidth);

legend("Clean", "Takeoff", "Landing", "Location", "Best")

xticks(0:0.02:0.6);      % Major ticks every 0.02 on X-axis
yticks(0:0.2:2.0);       % Major ticks every 0.2 on Y-axis

% Enable grid only at major ticks
ax = gca;                % Get current axes
ax.GridLineStyle = '-';  % Solid line for major grid

% Set minor ticks without grid lines
ax.XMinorTick = 'on';           % Enable minor ticks on X-axis
ax.YMinorTick = 'on';           % Enable minor ticks on Y-axis
ax.MinorGridLineStyle = 'none'; % Turn off minor grid lines

% Define minor tick intervals
ax.XAxis.MinorTickValues = 0:0.005:0.6;  % Minor ticks every 0.005 on X-axis
ax.YAxis.MinorTickValues = 0:0.05:2;   % Minor ticks every 0.05 on Y-axis

% Label axes
xlabel('C_D');           % Label for X-axis
ylabel('C_L', Rotation= 0);           % Label for Y-axis

% Set font size and other formatting adjustments to match style
ax.FontSize = fontSize_axes;        % Adjust font size
ax.LineWidth = 1;        % Set axis line width

ylim([0 2]);

% for shape of drag polar confirmation check Takahashi p. 135


%% functions for interpolation

function [lambda2] = lambda2fun(deflection, t_c, type)
%
if type == 'p' || type == "sp"
    % Import data sets
    data = importfileGeneral("lambda2.csv",4); % Call function to get first data set


    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2
    y1_3 = data(:, 4); % Fourth column is y data for curve 3

    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, deflection, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, deflection, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, deflection, 'linear', 'extrap');

    % Linearly interpolate between the three curves based on t_c

    if type == "sp"
        if t_c >=0.12 && t_c<=0.21
            % Interpolate between curve 1 and curve 2
            lambda2 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
        elseif t_c<0.12 || t_c >0.3
            error("t/c must be between 0.12 and 0.3");
        else
            % Interpolate between curve 2 and curve 3
            lambda2 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
        end
        % Output the interpolated results
    else
        lambda2 = y1_1_interp; %plain flap interpolation
    end
end

if type == "sl"
    % Import data sets
    data = importfileGeneral("lambda2_slotted.csv",4); % Call function to get first data set


    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % 0.21 Second column is y data for curve 1
    y1_2 = data(:, 3); % 0.30 Third column is y data for curve 2
    y1_3 = data(:, 4); % 0.12 Fourth column is y data for curve 3

    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, deflection, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, deflection, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, deflection, 'linear', 'extrap');

    % Linearly interpolate between the three curves based on t_c
    if deflection >=30
        if t_c >=0.21 && t_c<=0.30
            % Interpolate between curve 1 and curve 2
            lambda2 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
        elseif t_c<0.12 || t_c >0.3
            error("t/c must be between 0.12 and 0.3");
        else
            % Interpolate between curve 2 and curve 3
            lambda2 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
        end
        % Output the interpolated results

    else

        if t_c >=0.12 && t_c<=0.21
            % Interpolate between curve 1 and curve 3
            lambda2 = y1_1_interp + (y1_3_interp - y1_1_interp) * t_c / 0.5;
        elseif t_c<0.12 || t_c >0.3
            error("t/c must be between 0.12 and 0.3");
        else
            % Interpolate between curve 2 and curve 3
            lambda2 = y1_3_interp + (y1_2_interp - y1_3_interp) * (t_c - 0.5) / 0.5;
        end

    end

end

end

function [lambda3] = lambda3fun(bf_b)
%
    % Import data sets
    data = importfileGeneral("lambda3.csv",2); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
 
  
    % Perform interpolation for the first data set using C_LT as the query point
    lambda3 = interp1(x1, y1_1, bf_b, 'linear', 'extrap');
     
end


function [delta2] = delta2fun(deflection, t_c, type)
%
if type == "sp" || type == "p"
    % Import data sets
    data = importfileGeneral("delta2.csv",5); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2
    y1_3 = data(:, 4); % Fourth column is y data for curve 3
    y1_4 = data(:, 5); % Fifth column in y data for curve 4 (plain flap)
  
    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, deflection, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, deflection, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, deflection, 'linear', 'extrap');
    y1_4_interp = interp1(x1, y1_4, deflection, 'linear', 'extrap');

    if type == "sp"
        % Linearly interpolate between the three curves based on t_c

        if t_c >=0.21 && t_c<=0.30
            % Interpolate between curve 1 and curve 2
            delta2 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
        elseif t_c<0.12 || t_c >0.3
            error("t/c must be between 0.12 and 0.3");
        else
            % Interpolate between curve 2 and curve 3
            delta2 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
        end
        % Output the interpolated results
    end
    
    if type == 'p'
        delta2 = y1_4_interp;
    end
end

if type == "sl"
    % Import data sets
    data = importfileGeneral("delta2_slotted.csv",4); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2
    y1_3 = data(:, 4); % Fourth column is y data for curve 3
    
    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, deflection, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, deflection, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, deflection, 'linear', 'extrap');
    

    % Linearly interpolate between the three curves based on t_c

        if t_c >=0.12 && t_c<=0.21
            % Interpolate between curve 1 and curve 2
            delta2 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
        elseif t_c<0.12 || t_c >0.3
            error("t/c must be between 0.12 and 0.3");
        else
            % Interpolate between curve 2 and curve 3
            delta2 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
        end
        % Output the interpolated results

end
end

% function [delta3] = delta3fun(bf_b)
% %
%     % Import data sets
%     data = importfileGeneral("delta3.csv",2); % Call function to get first data set
% 
% 
%     % Extract x and y data from the first data set
%     % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
%     x1 = data(:, 1); % First column is x data
%     y1_1 = data(:, 2); % Second column is y data for curve 1
% 
% 
%     % Perform interpolation for the first data set using C_LT as the query point
%     delta3 = interp1(x1, y1_1, bf_b, 'linear', 'extrap');
% 
% end

function [mu2] = mu2fun(bf_b)
%
    % Import data sets
    data = importfileGeneral("mu2.csv",2); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
 
  
    % Perform interpolation for the first data set using C_LT as the query point
    mu2 = interp1(x1, y1_1, bf_b, 'linear', 'extrap');
     
end

function [K] = Kfun(bf_b)
%
    % Import data sets
    data = importfileGeneral("K.csv",2); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
 
  
    % Perform interpolation for the first data set using C_LT as the query point
    K = interp1(x1, y1_1, bf_b, 'linear', 'extrap');
     
end

function [FA_F6] = FA_F6fun(A)
%
    % Import data sets
    data = importfileGeneral("FA_F6.csv",2); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
 
  
    % Perform interpolation for the first data set using C_LT as the query point
    FA_F6 = interp1(x1, y1_1, A, 'linear', 'extrap');
     
end

function cf_c = cf_cfun(lambda1)
    % Import data sets
    data = importfileGeneral("lambda1.csv",2); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 2); % First column is x data
    y1_1 = data(:, 1); % Second column is y data for curve 1
 
  
    % Perform interpolation for the first data set using C_LT as the query point
    cf_c = interp1(x1, y1_1, lambda1, 'linear', 'extrap');

end

function [delta1] = delta1fun(cf_c, t_c, type)

if type == "p" || type == "sp"
%
    % Import data sets
    data = importfileGeneral("delta1.csv",4); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2
    y1_3 = data(:, 4); % Fourth column is y data for curve 3
  
    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, cf_c, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, cf_c, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, cf_c, 'linear', 'extrap');

    % Linearly interpolate between the three curves based on t_c

    if t_c >=0.12 && t_c<=0.21
        % Interpolate between curve 1 and curve 2
        delta1 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
    elseif t_c<0.12 || t_c >0.3
        error("t/c must be between 0.12 and 0.3");
    else
        % Interpolate between curve 2 and curve 3
        delta1 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
    end
    % Output the interpolated results

end

if type == "sl"
    % Import data sets
    data = importfileGeneral("delta1_slotted.csv",3); % Call function to get first data set
    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2

 % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, cf_c, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, cf_c, 'linear', 'extrap');


     if t_c >=0.12 && t_c<=0.30
        % Interpolate between curve 1 and curve 2
        delta1 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
     else
        error("thickness ratio must be between 0.12 and 0.30");
     end
end
end
function [delta3] = delta3fun(bf_b, lambda)
%
    % Import data sets
    data = importfileGeneral("delta3.csv",5); % Call function to get first data set
       

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2
    y1_3 = data(:, 4); % Fourth column is y data for curve 3
    y1_4 = data(:, 5); % Fifth column is y data for curve 4
    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, bf_b, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, bf_b, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, bf_b, 'linear', 'extrap');
    y1_4_interp = interp1(x1, y1_4, bf_b, 'linear', 'extrap');
    % Linearly interpolate between the three curves based on t_c

    if lambda <=1 && lambda>=0.5
        % Interpolate between curve 1 and curve 2
        delta3 = y1_1_interp + (y1_2_interp - y1_1_interp) * lambda / 0.5;
    elseif lambda<=0.5 && lambda >=0.33
        % Interpolate between curve 2 and curve 3
        delta3 = y1_2_interp + (y1_3_interp - y1_2_interp) * (lambda - 0.5) / 0.5;
    elseif lambda<= 0.33 && lambda >= 0.25
        % Interpolate between curve 3 and curve 4
        delta3 = y1_3_interp + (y1_4_interp - y1_3_interp) * (lambda - 0.5) / 0.5;
    else
        error('bf_b must be between 0.25 and 1');
    end
    % Output the interpolated results
end

function [mu1] = mu1fun(cf_c, t_c)
%
    % Import data sets
    data = importfileGeneral("delta3.csv",5); % Call function to get first data set
    

    % Extract x and y data from the first data set
    % Assume data has four columns: x, y for curve 1, y for curve 2, y for curve 3
    x1 = data(:, 1); % First column is x data
    y1_1 = data(:, 2); % Second column is y data for curve 1
    y1_2 = data(:, 3); % Third column is y data for curve 2
    y1_3 = data(:, 4); % Fourth column is y data for curve 3
    % Perform interpolation for the first data set using C_LT as the query point
    y1_1_interp = interp1(x1, y1_1, cf_c, 'linear', 'extrap');
    y1_2_interp = interp1(x1, y1_2, cf_c, 'linear', 'extrap');
    y1_3_interp = interp1(x1, y1_3, cf_c, 'linear', 'extrap');
    % Linearly interpolate between the three curves based on t_c

    if t_c <=0.3 && t_c>=0.21
        % Interpolate between curve 1 and curve 2
        mu1 = y1_1_interp + (y1_2_interp - y1_1_interp) * t_c / 0.5;
    elseif t_c<=0.21 && t_c >=0.12
        % Interpolate between curve 2 and curve 3
        mu1 = y1_2_interp + (y1_3_interp - y1_2_interp) * (t_c - 0.5) / 0.5;
    else
        error('bf_b must be between 0.12 and 0.3');
    end
    % Output the interpolated results
end