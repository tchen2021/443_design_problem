function output = importfileGeneral(filename, columns, dataLines)
%IMPORTFILE Import data from a text file
%  output = IMPORTFILETWOCOLUMN(FILENAME) reads data from text file FILENAME for the
%  default selection. For two column data Returns the numeric data.
%
%  output = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-columns array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  output = importfile("/Users/tysonc/Documents/GitHub/433_design_problem/drag_buildup/derivatives/B53.csv", 2 [2, Inf]);
%  This imports the data from the specified filepath, with two columns of
%  data, and from row 2 to the end (row 1 is the column labels x and y...)
%  See also READTABLE.


%% Input handling

% If dataLines is not specified, define defaults
if nargin < 3
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data


opts = delimitedTextImportOptions("NumVariables", columns);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

opts.VariableNames = ["x"];
opts.VariableTypes = ["double"];
VariableTypeAppend = "double";

 for i = 2:columns %need at least 2 columns (x,y pair)  
    VariableNameAppend = "VarName" + num2str(i);   
    % Specify column names and types
    opts.VariableNames = [opts.VariableNames(:, VariableNameAppend)];
    opts.VariableTypes = [opts.VariableTypes(:), VariableTypeAppend];
end

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
output = readtable(filename, opts);
%% Convert to output type
output = table2array(output);
end