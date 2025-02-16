% function output = importfileGeneral(filename, columns, dataLines)
% %IMPORTFILE Import data from a text file
% %  output = IMPMORTFILEGENERAL(FILENAME) reads data from text file FILENAME for the
% %  default selection. For two column data Returns the numeric data.
% %
% %  output = IMPORTFILE(FILE, DATALINES) reads data for the specified row
% %  interval(s) of text file FILENAME. Specify DATALINES as a positive
% %  scalar integer or a N-by-columns array of positive scalar integers for
% %  dis-contiguous row intervals.
% %
% %  Example:
% %  output = importfile("/Users/tysonc/Documents/GitHub/433_design_problem/drag_buildup/derivatives/B53.csv", 2, [2, Inf]);
% %  This imports the data from the specified filepath, with two columns of
% %  data, and from row 2 to the end (row 1 is the column labels x and y...)
% %  See also READTABLE.
% 
% 
% %% Input handling
% 
% % If dataLines is not specified, define defaults
% if nargin < 3
%     dataLines = [2, Inf];
% end
% 
% %% Set up the Import Options and import the data
% 
% 
% opts = delimitedTextImportOptions("NumVariables", columns);
% 
% % Specify range and delimiter
% opts.DataLines = dataLines;
% opts.Delimiter = ",";
% 
% opts.VariableNames = ["x"];
% opts.VariableTypes = ["double"];
% VariableTypeAppend = "double";
% 
% % for i = 2:columns %need at least 2 columns (x,y pair)  
% %     VariableNameAppend = "VarName" + num2str(i);   
% %     % Specify column names and types
% %     opts.VariableNames = [opts.VariableNames, VariableNameAppend];
% %     opts.VariableTypes = [opts.VariableTypes, VariableTypeAppend];;
% % end
% 
% for i = 2:columns  % Add column names and types dynamically
%     opts.VariableNames{end+1} = "VarName" + i;
%     opts.VariableTypes{end+1} = "double";
% end
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Import the data
% output = readtable(filename, opts);
% %% Convert to output type
% output = table2array(output);
% end

function output = importfileGeneral(filename, columns, dataLines)
%IMPORTFILEGENERAL Import data from a text file
%  output = IMPORTFILEGENERAL(FILENAME, COLUMNS) reads data from text file
%  FILENAME with the specified number of COLUMNS. The default row range
%  is from the second row to the end.
%
%  output = IMPORTFILEGENERAL(FILENAME, COLUMNS, DATALINES) reads data
%  for the specified row interval(s) of the text file. DATALINES should be
%  a scalar or an N-by-2 array specifying row intervals.
%
%  Example:
%  output = importfileGeneral("data.csv", 2, [2, Inf]);
%  This imports data from "data.csv" with two columns, starting from row 2.
%
%  See also: READTABLE, TABLE2ARRAY.

%% Input handling
if nargin < 3
    dataLines = [2, Inf]; % Default to reading from row 2 to the end
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", columns);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Initialize column names and types
opts.VariableNames = {'x'};  % Cell array of character vectors
opts.VariableTypes = {'double'};

for i = 2:columns  % Add column names and types dynamically
    opts.VariableNames{end+1} = ['VarName', num2str(i)];  % Convert to char
   
end
opts.VariableTypes = repmat({'double'}, 1, columns);

% Specify file-level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
output = readtable(filename, opts);

%% Convert to output type
output = table2array(output);
end