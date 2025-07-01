clear all 
close all
BST = load('C:\SortingFolder\HIP\October\spiketableKS4_MUA_october');
x1 = 250; %Start of window of analysis
x2 = 650; %End of window of analysis
%% 
for i = 1:size(BST.spiketable,1)
    % for k = 1:size(spiketable.oddData{1},1)
        [HistogramsSpikes] = SignificanceWindow_Spikes(BST.spiketable.oddData{i,1},x1,x2);
        
        HistogramsSpikesTable_MUA(i,:) = HistogramsSpikes;     
end


%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'C:\SortingFolder\HIP\October';
%%%%%%%%%%%%%%%%%%%%%%%%%

cd (savingPath);
save('HistogramsSpikesTable_MUA.mat', 'HistogramsSpikesTable_MUA');
%% 
