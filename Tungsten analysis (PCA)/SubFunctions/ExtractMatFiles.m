function [ BST, SPK, SPK9 ] = ExtractMatFiles(UnitTable,pathsSpikesPCA)
%Read and save in a structure the complete spikedata, spikedata9 and the
%hole BST
% spikedata9 just includes trials that have 9 standard before the deviant
% in the oddball paradigm, and in the case of the controls ans periodic
% oddball it takes the first 9 sequences. 

for j = 1:size(UnitTable,1)
    cd (pathsSpikesPCA(j,1))
    Files = dir(fullfile(pathsSpikesPCA(j,1), '*.mat'));
    fieldname = UnitTable(j,6).Name;

    for k = 1:size(Files,1)

        if contains(Files(k).name,'bst') && contains(Files(k).name,UnitTable(j,5).Spike)
            PathBST = join([pathsSpikesPCA(j,1),"\",Files(k).name]);
            PathBST = erase(PathBST, " ");
            tmpBST = load(PathBST);
            BST(j,1)= tmpBST;
        end
        
          if contains(Files(k).name,'spikedata') && contains(Files(k).name,UnitTable(j,5).Spike) && strlength((Files(k).name))<=15
            PathSPK = join([pathsSpikesPCA(j,1),"\",Files(k).name]);
            PathSPK = erase(PathSPK, " ");
            tmpSPK = load(PathSPK);
            SPK(j,1) = tmpSPK;
          end

         if contains(Files(k).name,'spikedata') && contains(Files(k).name,'_9') && contains(Files(k).name,UnitTable(j,5).Spike)
            PathSPK9 = join([pathsSpikesPCA(j,1),"\",Files(k).name]);
            PathSPK9 = erase(PathSPK9, " ");
            tmpSPK9 = load(PathSPK9);
            SPK9(j,1)= tmpSPK9;

         end
    end
end

if size(UnitTable,1) ~= size(SPK9,1)
    error('There are less files than neurons in the Unit Table')
end

end