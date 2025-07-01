function [ BST, BST_MUA ] = ExtractDataMultiChannel(UnitTable,pathsSpikes)
%Read and save in a structure the complete spikedata, spikedata9 and the
%hole BST
% spikedata9 just includes trials that have 9 standard before the deviant
% in the oddball paradigm, and in the case of the controls ans periodic
% oddball it takes the first 9 sequences. 

for j = 1:size(UnitTable,1)
    cd (pathsSpikes(j,1))
    Files = dir(fullfile(pathsSpikes(j,1), '*.mat'));
    fieldname = UnitTable(j,1).Animal;

    for k = 1:size(Files,1)

        if contains(Files(k).name,'bst_KS4') && length(Files(k).name) == 11
            PathBST = join([pathsSpikes(j,1),"\",Files(k).name]);
            PathBST = erase(PathBST, " ");
            tmpBST = load(PathBST);
            BST(j,1)= tmpBST;
        end
        
          if contains(Files(k).name,'bst_KS4_mua')
            PathBST_MUA = join([pathsSpikes(j,1),"\",Files(k).name]);
            PathBST_MUA = erase(PathBST_MUA, " ");
            tmpBST_MUA = load(PathBST_MUA);
            BST_MUA(j,1) = tmpBST_MUA;
          end

    end
end

if size(UnitTable,1) ~= size(BST,1)
    error('There are less files than neurons in the Unit Table')
end

end