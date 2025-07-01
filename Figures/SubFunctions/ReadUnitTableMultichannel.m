function [UnitTable] = ReadUnitTableMultichannel(Rango,Page,Path)
%Reads an excel spreadsheet with the data of unit table
opts = spreadsheetImportOptions("NumVariables", 18);
% Specify sheet and range
opts.DataRange = Rango;
opts.Sheet = Page;
% Specify column names and types
opts.VariableNames = ["Animal", "Tank", "Tract", "Depth", "Coordinates"];
opts.VariableTypes = ["string","string","double","double","string"];

% Import table
UnitTable  = readtable(Path, opts);

clear opts

end