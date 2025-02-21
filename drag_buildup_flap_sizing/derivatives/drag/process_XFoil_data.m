function [value1, value2] = process_XFoil_data(airfoil, lookup_value, Re_target, lookup_type)
    addpath XFoil_Results
    % Folder where the XFoil results are stored
    outputFolder = 'XFoil_Results';
    
    % Get list of existing airfoil files
    airfoilFiles = dir(fullfile(outputFolder, ['XFoil_' airfoil '*.mat']));
    
    if isempty(airfoilFiles)
        error('Airfoil file does not exist.');
    end

    % Extract Reynolds numbers from the file names
    Re_values = [];
    for i = 1:length(airfoilFiles)
        % Use regular expression to extract Reynolds number from filename
        reMatch = regexp(airfoilFiles(i).name, 'Re(\d+)', 'tokens');
        if ~isempty(reMatch)
            Re_values = [Re_values, str2double(reMatch{1}{1})];
        end
    end
    
    % Remove duplicates and sort the Reynolds numbers
    Re_values = unique(Re_values);
    
    % Display available Reynolds numbers
%     disp('Available Reynolds numbers:');
%     disp(Re_values);
    
    % Initialize output variables
    value1 = NaN;
    value2 = NaN;

    % Check if the specified Reynolds number is available
    if ismember(Re_target, Re_values)
        % Reynolds number exists, load the file and return data
        fileName = airfoilFiles(find(Re_values == Re_target, 1)).name;
%         disp(['Loading file: ', fileName]);
        load(fullfile(outputFolder, fileName), 'pol', 'foil');
        
        % Extract the data
        CL = pol.CL;
        CD = pol.CD;
        alpha = pol.alpha;
        
        % Determine the lookup type
        if strcmpi(lookup_type, 'CL')
            [~, idx] = min(abs(CL - lookup_value));
            value1 = CD(idx);
            value2 = alpha(idx);
%             disp(['Closest CL: ', num2str(CL(idx)), ', Closest CD: ', num2str(CD(idx)), ', Closest alpha: ', num2str(alpha(idx))]);
        elseif strcmpi(lookup_type, 'alpha')
            [~, idx] = min(abs(alpha - lookup_value));
            value1 = CL(idx);
            value2 = CD(idx);
%             disp(['Closest alpha: ', num2str(alpha(idx)), ', Closest CL: ', num2str(CL(idx)), ', Closest CD: ', num2str(CD(idx))]);
        else
            error('Invalid lookup type. Use ''CL'' or ''alpha''.');
        end
        
    else
        % Reynolds number does not exist, interpolate between the closest values
        Re_values = sort(Re_values);
        lowerRe = max(Re_values(Re_values <= Re_target));
        upperRe = min(Re_values(Re_values >= Re_target));
        
        % If target Reynolds number is outside the available range
        if isempty(lowerRe) || isempty(upperRe)
            error('Reynolds number out of range.');
        end
        
        % Load data for both Reynolds numbers
        lowerFile = airfoilFiles(find(Re_values == lowerRe, 1)).name;
        upperFile = airfoilFiles(find(Re_values == upperRe, 1)).name;
        
%         disp(['Loading lower file: ', lowerFile]);
        load(fullfile(outputFolder, lowerFile), 'pol', 'foil');
        lowerData = struct('CL', pol.CL, 'CD', pol.CD, 'alpha', pol.alpha);
        
%         disp(['Loading upper file: ', upperFile]);
        load(fullfile(outputFolder, upperFile), 'pol', 'foil');
        upperData = struct('CL', pol.CL, 'CD', pol.CD, 'alpha', pol.alpha);
        
        % Interpolate based on lookup type
        if strcmpi(lookup_type, 'CL')
            [~, lowerIdx] = min(abs(lowerData.CL - lookup_value));
            [~, upperIdx] = min(abs(upperData.CL - lookup_value));
            lowerValue1 = lowerData.CD(lowerIdx);
            upperValue1 = upperData.CD(upperIdx);
            lowerValue2 = lowerData.alpha(lowerIdx);
            upperValue2 = upperData.alpha(upperIdx);
        elseif strcmpi(lookup_type, 'alpha')
            [~, lowerIdx] = min(abs(lowerData.alpha - lookup_value));
            [~, upperIdx] = min(abs(upperData.alpha - lookup_value));
            lowerValue1 = lowerData.CL(lowerIdx);
            upperValue1 = upperData.CL(upperIdx);
            lowerValue2 = lowerData.CD(lowerIdx);
            upperValue2 = upperData.CD(upperIdx);
        else
            error('Invalid lookup type. Use ''CL'' or ''alpha''.');
        end
        
        % Interpolate the values
        value1 = interp1([lowerRe, upperRe], [lowerValue1, upperValue1], Re_target, 'linear');
        value2 = interp1([lowerRe, upperRe], [lowerValue2, upperValue2], Re_target, 'linear');
        
%         % Display interpolated values
%         if strcmpi(lookup_type, 'CL')
%             disp(['Interpolated CL: ', num2str(lookup_value), ', Interpolated CD: ', num2str(value1), ', Interpolated alpha: ', num2str(value2)]);
%         elseif strcmpi(lookup_type, 'alpha')
%             disp(['Interpolated alpha: ', num2str(lookup_value), ', Interpolated CL: ', num2str(value1), ', Interpolated CD: ', num2str(value2)]);
%         end
    end
end
