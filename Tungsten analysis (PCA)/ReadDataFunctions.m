clear all

[pathsSpikes tmp] = ReadPCAPaths('D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\HIPDATA\spikes');

[HIP] = ReadUnitTable("A2",3,'D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\HIPDATA\Unit_Table_HIP_Simplified_Enero_2024.xlsx');

[pathsSpikesPCA] = pathsPCA(tmp,pathsSpikes,HIP);

[BST, SPK, SPK9] = ExtractMatFiles(HIP,pathsSpikesPCA);
y1 = 250;
y2 = 650;
%%

for i = 1:size(HIP,1)
[NormalizedIndexes] = IndexGeneratorDensity(SPK(i).spikedata,HIP(i,7).x1,HIP(i,8).x2);
NormalizedIndexesTable(i,:) = NormalizedIndexes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
savingPath = 'D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\HIPDATA\Revised\Sub';
%%%%%%%%%%%%%%%%%%%%%%%%%

cd (savingPath);
save('NormalizedIndexesTable.mat', 'NormalizedIndexesTable');

%%
for i = 1:size(HIP,1)
[NormalizedIndexesC] = IndexGeneratorDensity(SPK(i).spikedata,y1,y2);
NormalizedIndexesTableC(i,:) = NormalizedIndexesC;
end

cd (savingPath);
save('NormalizedIndexesTableC.mat', ['NormalizedIndexesTableC']);

%%
for i = 1:size(HIP,1)
[NormalizedIndexes9] = IndexGeneratorDensity(SPK9(i).spikedata9,HIP(i,7).x1,HIP(i,8).x2);
NormalizedIndexesTable9(i,:) = NormalizedIndexes9;
end

cd (savingPath);
save('NormalizedIndexesTable9.mat', ['NormalizedIndexesTable9']);

%%
for i = 1:size(HIP,1)
[NormalizedIndexes9C] = IndexGeneratorDensity(SPK9(i).spikedata9,y1,y2);
NormalizedIndexesTable9C(i,:) = NormalizedIndexes9C;
end

cd (savingPath);
save('NormalizedIndexesTable9C.mat', ['NormalizedIndexesTable9C']);