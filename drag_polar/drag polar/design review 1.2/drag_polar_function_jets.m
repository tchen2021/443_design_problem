function [Range_LD_0, Range_LD_1, Endurance_LD_0, Endurance_LD_1, CL_Endurance_0, CL_Endurance_1] = drag_polar_function_jets
%Function-format of drag_polar script by Tyson Chen
%requires dragpolar_data.mat in the same folder

%OUTPUTS:
%Range_LD_0/1: max range L/D at drag index 0/1
%Endurance_LD_0/1: max endurance L/D at drag index 0/1
%CL_Endurance_0/1: CL at endurance condition at index 0/1

%% load data and L/D from competitive asssessment
dragpolar_data = load('dragpolar_data.mat');
CL_market = dragpolar_data.ComparisionTablenov4S2.CL; %CL data
CD_market = dragpolar_data.ComparisionTablenov4S2.CD; %CD data
Market_names = dragpolar_data.ComparisionTablenov4S2.Name;
% i think the range should be about 8 to 13 - from sam 

%% first drag polar (props)
% Example point 
%CL_dat = 0.40;
%CD_dat = 0.025 + 0.01;


LD(1) = 12.15224; % T-6B 

x(:,1) = linspace(0, 0.1);
y(:,1) = LD(1).*x(:,1); % this is just y = mx + b; to plot the line 


% recall the equation for a drag polar 
% CDtot = CD0 + 1/k * CL^2 

% pick your range of CL, at this point its just a guess 
CL(:,1) = linspace(0, 1.5); 
CD0(1) = 0.015; % also just a guess 

e(1) = 0.70; % guess with this one too to get it close  %0.3
A(1) = 6.515662651; %T-6B 
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

LD(2) = 9.722; % OV-10D 

x(:,2) = linspace(0, 0.15);
y(:,2) = LD(2).*x(:,2); % this is just y = mx + b; to plot the line 


% recall the equation for a drag polar 
% CDtot = CD0 + 1/k * CL^2 

% pick your range of CL, at this point its just a guess 
CL(:,2) = linspace(0, 1.5); 
CD0(2) = 0.035; % also just a guess 

e(2) = 0.7; % guess with this one too to get it close  0.21
A(2) = 2.8; %OV-10 (defunct) 
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
%% drag polars (jets)

LD(3) = 13.3675; % A-37 Cessna 

x(:,3) = linspace(0, 0.15);
y(:,3) = LD(3).*x(:,3); % this is just y = mx + b; to plot the line 

LD(4) = 11.6; % L-39 Albotross 

x(:,4) = linspace(0, 0.15);
y(:,4) = LD(4).*x(:,4); % this is just y = mx + b; to plot the line 


% recall the equation for a drag polar 
% CDtot = CD0 + 1/k * CL^2 

% pick your range of CL, at this point its just a guess 
CL(:,3) = linspace(0, 1.5); 
CD0(3) = 0.02; % also just a guess 

e(3) = 0.17; % guess with this one too to get it close 
A(3) = 6.2; 
k(3) = (pi * e(3) * A(3)); 


CDtot(:,3) = CD0(3) + 1/k(3) * CL(:,3).^2; 

CL(:,4) = linspace(0, 1.5); 
CD0(4) = 0.08; % also just a guess 

e(4) = 0.055; % guess with this one too to get it close 
A(4) = 4.8; %
k(4) = (pi * e(4) * A(4)); 


CDtot(:,3) = CD0(3) + 1/k(3) * CL(:,3).^2; 

CDtot(:,4) = CD0(4) + 1/k(4) * CL(:,4).^2; 
%gradient2 = gradient(CL(:,2)/CDtot(:,2));

%find the max L/D CL point, also the max endurance 
[placeholder, idx_max(3)] = max(CL(:,3)./CDtot(:,3))
[placeholder, idx_max(4)] = max(CL(:,4)./CDtot(:,4))
CL_max(3) = CL(idx_max(3),3);
CD_max(3) = CDtot(idx_max(3),3);
LD_max(3) = CL_max(3)/CD_max(3);
CL_max(4) = CL(idx_max(4),4);
CD_max(4) = CDtot(idx_max(4),4);
LD_max(4) = CL_max(4)/CD_max(4);

%find the max range for jet
[placeholder, idx_range(3)] = max((CL(:,3)).^0.5 ./ CDtot(:,3))
[placeholder, idx_range(4)] = max((CL(:,4)).^0.5 ./ CDtot(:,4))
CL_range(3) = CL(idx_range(3),3);
CD_range(3) = CDtot(idx_range(3),3);
LD_range(3) = CL_range(3)/CD_range(3);
CL_range(4) = CL(idx_range(4),4);
CD_range(4) = CDtot(idx_range(4),4);
LD_range(4) = CL_range(4)/CD_range(4);
%find the max endurance for jet
[placeholder, idx_max_endurance(3)] = max((CL(:,3)) ./ CDtot(:,3))
[placeholder, idx_max_endurance(4)] = max((CL(:,4)) ./ CDtot(:,4))
CL_endurance(3) = CL(idx_max_endurance(3),3);
CD_endurance(3) = CDtot(idx_max_endurance(3),3);
LD_endurance(3) = CL_endurance(3)/CD_endurance(3);
CL_endurance(4) = CL(idx_max_endurance(4),4);
CD_endurance(4) = CDtot(idx_max_endurance(4),4);
LD_endurance(4) = CL_endurance(4)/CD_endurance(4);


Range_LD_0 = LD_range(3);
Range_LD_1 = LD_range(4); 
Endurance_LD_0 = LD_endurance(3);
Endurance_LD_1 = LD_endurance(4);

CL_Endurance_0 = CL(idx_max_endurance(3),3);
CL_Endurance_1 = CL(idx_max_endurance(4),4);

end

