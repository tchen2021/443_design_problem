function [CD0W] = iterator(airfoil,alpha,iprime_r, epsilon,cbarbari, S_arr,V,rho,Se, S,h,series)
%ITERATOR Summary of this function goes here
%creates an array with rows of 

%   Detailed explanation goes here
% airfoil = 'NACA 0012';
% lookup_type = 'alpha'; % Options: 'CL' or 'alpha'
% lookup_value = 5.0; % Target value for CL or alpha
% Re_target = 2e6;


%check inputs are the same size
% if length(iprime_r) ~= length(S_arr) || length(epsilon) ~= length(S_arr)
% 
%     error('dimensions of either alpha, incidence angle, or wing section areas not equal in length!');
% end
V = V * 1.68781;        %convert knots to ft/s
CD0W_arr = zeros(length(S_arr));

mu = SutherlandsEquation(h);

lookup_type = 'alpha'; % Options: 'CL' or 'alpha'
epsilon=epsilon*ones(1,length(S_arr));
iprime_r=iprime_r*ones(1,length(S_arr));


if strcmpi(lookup_type, 'alpha')
    
    for i=1:length(S_arr)
    Re_i = (V* rho * cbarbari(i) )/ mu;
    disp(Re_i);
    alpha_i = alpha + iprime_r(i) + epsilon(i);     %should be in degrees
    [~, CD0W_arr(i)] = process_XFoil_data(airfoil, alpha_i, Re_i, 'alpha', series);
    %disp(['Closest CL: ', num2str(CD0W(i)), ', Closest CD: ', num2str(result2)]);
    end
 
    CD0W = 2 * (Se/S) * sum(sum(CD0W_arr .* S_arr') / sum(S_arr));
else
    error('Invalid lookup type. Use ''CL'' or ''alpha''.');
end


end

