function [UnitTable] = ReadUnitTable(Rango,Page,Path)
%Reads an excel spreadsheet with the data of unit table
opts = spreadsheetImportOptions("NumVariables", 18);
% Specify sheet and range
opts.DataRange = Rango;
opts.Sheet = Page;
% Specify column names and types
opts.VariableNames = ["Animal", "Unit", "Division", "ThreshSTDs", "Spike","Name","x1","x2"];
opts.VariableTypes = ["string","string","string","double","string","string","double","double"];

% Import table
UnitTable  = readtable(Path, opts);

clear opts

end