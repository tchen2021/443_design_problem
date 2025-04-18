function dragdelta2 = importfiledragdelta2(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DRAGDELTA2 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  DRAGDELTA2 = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  dragdelta2 = importfile("/Users/tysonc/Documents/GitHub/433_design_problem/drag_buildup/derivatives/drag_delta2.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 31-Jan-2025 01:55:51

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["x", "Curve1"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
dragdelta2 = readtable(filename, opts);

%% Convert to output type
dragdelta2 = table2array(dragdelta2);
end