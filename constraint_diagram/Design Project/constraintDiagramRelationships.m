%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint Diagram %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pick which loading Method you desire
    % load planeData
    fileLocation="C:\Users\bravo_4e3\OneDrive - Cal Poly\AERO 435 - Design I\Design Project\planeData.xlsx";
    planeData = importplaneData(fileLocation, "Sheet1", [2, Inf])
    
    % Classify Data8uj 
    names=planeData{:,1};
    MTOW=planeData{:,2};
    W_e=planeData{:,3};
    S=planeData{:,4};
    AR=planeData{:,5};
    V_s=planeData{:,6};
    V_max=planeData{:,7};
    P=planeData{:,8};
    S_TOFL=planeData{:,9};
    S_L=planeData{:,10};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumptions  Overview
h=0;% altitude range [ft]
% Assume Cl_MaxTO = Cl_MaxL


%%% Standard Atmo Values %%%
[p0, T0, rho0, a0] = atmosphere(0); % Calculate standard atmosphere values at sea level
[p, T, rho, a] = atmosphere(h);% Calculate standard atmosphere values at variable height
% Correct for non standard atmosphere
sigma=rho./rho0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Graph STOFL vs TOP
for i =1: height(planeData)
    % Quick Conversions
    W=MTOW(i);
    WS=W/S(i);

    %% Calculate Takeoff Parameters & Graph - Got rid of assumed ratio
    %%% Assumption - Assume Cl_MaxTO = Cl_MaxL
    % Calculate Cl_maxL from stall speed
    Cl_maxL=(2*W)/(rho*S(i)*V_s(i)^2);
    % CL_max TO 
    Cl_maxTO=Cl_maxL;
    %%% Calculate TOP (Takeoff Parameters) %%%
    TOP(i)=1/(WS / ( (P(i)/W) *Cl_maxTO*sigma));
    
    plot(1./TOP(i),S_TOFL(i),'*'); title("STOFL vs TOP - find best line"); hold on;

end
 grid on; xlabel("TOP"); ylabel("STOFL (ft)"); legend(names)


%% Plot the Landing distance
for i =1: height(planeData)
    figure (2)
    plot(V_s(i),S_L(i),"*"); title("Landing distance vs Stall speed"); hold on;
end
grid on; xlabel("Stall Speed (ft/s)"); ylabel("S_L (ft)"); legend(names)

