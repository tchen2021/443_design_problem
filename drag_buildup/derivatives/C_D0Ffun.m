function C_D0F = C_D0Ffun(Re, S, S_A_wet, S_B_wet, S_C_wet, S_CAB, flow_parameter, l_f, D, Lambda_f)
%assumes aerodynamically smooth surface finishing

E_e = l_f/D;        %slenderness of the fuselage [dimensionless]

%check if flow parameter is inputted correct
if flow_parameter == 't'
    data1 = importfiledrag_turbulent("drag_turbulent.csv");
    x1 = data1(:, 1); % First column is x data
    y1 = data1(:, 2); % Second column is y data for curve 1
    C_f = interp1(x1, y1, Re, 'linear', 'extrap');
elseif flow_parameter == 'l'
    data2 = importfiledrag_turbulent("drag_laminar.csv");
    x2 = data2(:, 1); % First column is x data
    y2 = data2(:, 2); % Second column is y data for curve 1
    C_f = interp1(x2, y2, Re, 'linear', 'extrap');
else
    error("Input must be 't' for tubulent and 'l' for laminar");
end


data3 = importfiledrag_turbulent("drag_F.csv");
x3 = data3(:, 1); % First column is x data
y3 = data3(:, 2); % Second column is y data for curve 1
F = interp1(x3, y3, E_e, 'linear', 'extrap');

data4 = importfiledrag_F("drag_K.csv");
x4 = data4(:, 1); % First column is x data
y4 = data4(:, 2); % Second column is y data for curve 1
K = interp1(x4, y4, Lambda_f, 'linear', 'extrap');

C_D0A = C_f * F * S_A_wet / S;                          %drag in region A
C_D0B = C_f * S_B_wet / S;                              %drag in region B
C_D0C = C_f * F *S_C_wet / S;                           %drag in region C
deltaC_D0LambdaF = K * C_D0C / 100;                     %additional drag due to tail up-sweep angle of the fuselage
deltaC_DS = 0.002;                                      %cabin/cockpit drag. Assumed to be curved windscreen with round upper edge. From Gudmundsson p.727
deltaC_D0CAB = deltaC_DS * (S_CAB/S);                   %additional drag due to cabin/cockpit protuburance.
C_D0F = C_D0A + C_D0B + C_D0C + deltaC_D0LambdaF + deltaC_D0CAB; %fuselage parasitic drag


end