close all 
clear 
clc 

addpath(genpath('C:\Users\diego\MATLAB\Tools'))
load("colormap_tab10.mat")
fontSize = 12; 

% alpha_f = Beta - sigma; % where sigma is sidewash 
d_sigma_d_Beta = 0; % sidewash derivative is 0 

%% Aircraft Parameters 
Sref = 275; 
Ss = 246; % Fuselage Side Reference Area (CHECK ME) 
b = 42.3; 

% Vertical Stab Parameters 
ARv = 1.8; 
bv = 7.2; 
cv = 4.33; 
Sv = 30; 

% Aircraft configuration 
idx = 2; % cr_c: [20 30 40] 
jdx = 3; % number of prop blades: [2 3 4 6] 

% Rudder Parameters 
cr_c = [0.20 0.30 0.40]; 
cr = cv * cr_c(idx); 
Sr = Sv*cr_c(idx); 


%% Set Environmental Conditions 
H = 0; 
[P, T, rho, a] = atmosphere(H);
%% Vertical Tail Contribution 
C_Vv = 0.07; 
VF_V = 1.1; 
af  = 2.8325; 

Cn_Betaf = C_Vv * af * (VF_V)^2 * (1 - d_sigma_d_Beta);  % f is "fin";

%% Wing Contribution 
Cn_BetaW = 0.02;  % F B 12,1 AR = 6.5 Lambda = 0 

%%  Fuselage Contribution
Lf = 34; 
d = 14.38; 
d_Lf = d/Lf;
Lf_b = Lf/b; 

DeltaCnBetaWB = 0; % low wing; 0.005 mid wing; 0.010 high wing

h1 = 4.53; 
w1 = 3.84; 
h2 = 3.9; 
w2 = 2.27; 

h = 4.2; 
Lf_h = Lf/h;
K_Beta = 0.125; 

Cn_BetaF = -0.96*K_Beta*(Ss/Sref) * (Lf_b) * (h1/h2)^(1/2) * (w1/w2)^(1/3);

%% Propulsion Contribution
lp = d - 2.3; 
N = 1; % number of props 
D = 7.5; % diameter (ft) 

dCyprop_dBeta = [-2.87e-5 -4.10e-5 -5.16e-5 -8.9e-5]; % number prop blades: [2 3 4 6] 
dCyprop_dBeta = dCyprop_dBeta(jdx);

Cn_Beta_prop = ((pi*D^2*lp)/(4*Sref*b)) * N * dCyprop_dBeta; % where N is number of propellers 





%% Directional Control 
ar = [1.3 1.6 1.8]; % 20 30 40 % of chord 
% Cn_delta_r = -ar(idx)*C_Vv*(VF_V)^2; % Rudder Power 
% Cn_delta_r = -ar(idx)*C_Vv*0.85; 
Cn_delta_r = -0.065; % - 0.005 - 0.02

%% Find Total Cn_Beta Derivative 

Cn_Beta = Cn_Betaf + Cn_BetaW + Cn_BetaF + Cn_Beta_prop  + Cn_delta_r; 

%% Rudder Force
b1 = [0.1998 0.244 0.2695]; 
b2 = [1.281 1.6554 1.985]; 

G = 0.099/1; % delta_e/delta_s

Pr_max_t = 150; % lbf 
Pr_max_s = 20; % lbf 

%% Define Velocity Limits
V_low = 110; 
V_high = 650; 

Vstall = 139; % ft/s 
Vmax = 506; % ft/s 

nfine = 100; 

fps2kts = 1/1.688; 

%% Fine Velocity Array 
% Define array of velocities
VF = linspace(V_low, V_high, nfine); 

% Preallocate arrays for results
Beta_Prmax_s = zeros(size(VF));
Beta_Prmax_t = zeros(size(VF));
Beta_Prmax_s_deg = zeros(size(VF));
Beta_Prmax_t_deg = zeros(size(VF));

% Sustained Limit 
for i = 1:length(VF)
    Beta_Prmax_s(i) = - 1 / (b1(idx) * (1 - d_sigma_d_Beta) + b2(idx) * (Cn_Beta / Cn_delta_r)) * ...
        (Pr_max_s / (G * (rho / 2) * VF(i)^2 * Sr * cr));
    Beta_Prmax_s_deg(i) = rad2deg(Beta_Prmax_s(i));
end

% Temporary Limit
for i = 1:length(VF)
    Beta_Prmax_t(i) = - 1 / (b1(idx) * (1 - d_sigma_d_Beta) + b2(idx) * (Cn_Beta / Cn_delta_r)) * ...
        (Pr_max_t / (G * (rho / 2) * VF(i)^2 * Sr * cr));
    Beta_Prmax_t_deg(i) = rad2deg(Beta_Prmax_t(i));
end


%% Directional Control 
% alpha_F = Beta - 0; 
% delta_r = -b1(idx)/b2(idx) * alpha_F;

delta_r = deg2rad(30);
Beta =  -Cn_delta_r/Cn_Beta * delta_r;

Beta_Limit_deg = rad2deg(Beta);


%% Velocity Conversion 
VF_kts = VF .* fps2kts; 
Vstall_kts = Vstall .* fps2kts; 
Vmax_kts = Vmax .* fps2kts; 

%% Figure 1 (BETA Limit) 

lineWidthVal = 1.7; 

f = figure ;
hold on 
grid on 
p1 = plot(VF_kts, Beta_Prmax_s_deg); 
p1.LineWidth = lineWidthVal; 
p1.Color = NewPolyGreen; 

p2 = plot(VF_kts, Beta_Prmax_t_deg); 
p2.LineWidth = lineWidthVal; 
p2.Color = Gold; 

p3 = xline(Vstall_kts); 
p3.LineWidth = lineWidthVal;

p4 = xline(Vmax_kts); 
p4.LineWidth = lineWidthVal;

p5 = plot(VF_kts, Beta_Limit_deg * ones(size(VF)));
p5.LineWidth = lineWidthVal;
p5.Color = BrickRed; 



% % Move the x-axis to the top
% ax = gca;
% ax.XAxisLocation = 'top';

ylim([0 50])
ylabel('\beta  [deg]')
xlabel('V [kts]')

set(gca,'FontSize', fontSize)

% set(f, 'Position', [100, 100, 1200, 600]);

%% Saving Figure 

% % Define the folder path and SVG file name
% folder_path = 'C:\Users\diego\MATLAB\Winter 2025\AERO 444\FIGURES';
% file_name = 'Beta_Limit.svg';
% 
% % Join the folder path and file name
% % Save the figure as an SVG
% saveas(gcf, fullfile(folder_path, file_name));

%% Crosswind 
Vv = zeros(size(VF));

% based on Max Rudder Def 
for i = 1:length(VF)
    Vv(i) = VF(i) * tan(Beta);
end

Vv_kts = Vv * fps2kts; 

% Sustained Limit 
Vv_Prmax_s = VF .* tand(Beta_Prmax_s_deg);
Vv_Prmax_s_kts = Vv_Prmax_s .* fps2kts; 

% Temporary Limit 
Vv_Prmax_t = VF .* tand(Beta_Prmax_t_deg);
Vv_Prmax_t_kts = Vv_Prmax_t .* fps2kts; 

%% Figure 2 (Vv Limit) 

f2 = figure;
hold on 
grid on 

p1 = plot(VF_kts(1:end), Vv_Prmax_s_kts(1:end)); 
p1.LineWidth = lineWidthVal; 
p1.Color = NewPolyGreen; 

p2 = plot(VF_kts(30:end), Vv_Prmax_t_kts(30:end)); 
p2.LineWidth = lineWidthVal; 
p2.Color = Gold; 

p3 = xline(Vstall_kts); 
p3.LineWidth = lineWidthVal;

p4 = xline(Vmax_kts); 
p4.LineWidth = lineWidthVal;

p5 = plot(VF_kts, Vv_kts);
p5.LineWidth = lineWidthVal;
p5.Color = BrickRed; 


ylabel('V_v [kts]')
xlabel('V [kts]')
ylim([0 125])

set(gca,'FontSize', fontSize)

% set(f2, 'Position', [100, 100, 1200, 600]);

%% Saving Figure 

% Define the folder path and SVG file name
file_name = 'CrossWind_Limit.svg';

% Join the folder path and file name
% Save the figure as an SVG
saveas(gcf, fullfile(folder_path, file_name));




 %% romans atmo
 function [p, T, rho, a] = atmosphere(h)
% atmosphere(h) function compute p, T, rho and a,  giving the altitude h 
% p : pressure in lb/ft^2
% T : temperature in degree Rankine
% rho : density in slug/ft^3
% a : speed of sound
% h : altitude in ft, h<=104990ft 

%% Constants
Cv = 4290;      % constant volume specific heat for the air
R = 1716;       % ideal gas constant
Cp = R+Cv;      % constant pressure specific heat for air
gamma = Cp/Cv;  % specific heat ratio
T_sl = 518.69;  % sea level temperature
p_sl = 2116.2;  % sea level pressure
%% temperature ratio
if  h<=36089
    theta = 1-(6.875*10^-6)*h;
elseif h>36089 && h<=65617
    theta = 0.75189;
elseif h>65617 && h<=104990
    theta = 0.75189+(1.0577*(10^-6))*(h-65617);
else
    theta=nan;
    warning('altitude exceeded')
end
        
%% pressure density

delta = ( 1-(6.875*(10^-6)).*h ).^5.2561.* ((0 <= h) & (h <= 36089)) +...
        (0.2234*exp(4.806*(10^-5).*(36089-h))).* ((36089 < h )& (h <= 65617)) +...
        (3.174716*(10^-6).*(0.75189+1.0577.*(10^-6).*(h-65617)).^-34.164).* ((65617<h) & (h<=104990));
%% p, T and rho

T = T_sl.*theta;
p = p_sl.*delta;
rho = p./(R.*T);
a = sqrt(gamma*R.*T);

 end




