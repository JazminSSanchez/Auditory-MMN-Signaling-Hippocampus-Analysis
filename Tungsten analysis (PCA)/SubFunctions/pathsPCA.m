function [pathsSpikesPCA] = pathsPCA(tmp,pathsSpikes,UnitTable)
%Creates a table with all the directories of interest from this area of the
%Hippocampus

for ii = 1:numel(tmp)
    for iii = 1:size(UnitTable,1)
        if tmp(ii).name == UnitTable.Animal(iii)
             pathsSpikesPCA(iii,1) = join([pathsSpikes(ii,1),"\",UnitTable.Unit(iii)]);
             pathsSpikesPCA(iii,1) = erase(pathsSpikesPCA(iii,1) ," ");
        end
    end
end

end