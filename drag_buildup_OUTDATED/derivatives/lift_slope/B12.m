function [CL_alpha] = B12(xaxis)
%Code figures out based on x value what the y value that matches with it
%should be

%% Fnd point on graph (B11B)
% load digitized graph data
B12 = importB12("B12-1.csv");
[~,Index] = min(abs(xaxis-B12(:,1))); % Find y from line according to x value 
CL_alpha=B12(Index,2);

end % end of function