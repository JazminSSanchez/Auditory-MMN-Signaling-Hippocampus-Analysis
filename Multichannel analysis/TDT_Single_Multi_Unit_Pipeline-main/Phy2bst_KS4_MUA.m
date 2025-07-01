%% Convert Sorted Phy data back to usable TDT (bst) data
% For TDT data that was prepared for Spike'Interface using StreamToBinary_AH.m
% Uses the sample numbers of where streams were appended to separate them
% again (needs the associated StreamSplitInfo.mat file).
%
% This converts Phy files from SpikeInterface to bst type data to match TDT epochs in BST files.
%
% Step 1 - import phy spike times and clusters and separate spikes
% by which TDT block they were part of.
% Step 2 - Convert Phy spike data to bst format, to make it comparable to
% TDT stimuli.
%
% AH 02/2023
% JSS 01/2025

clear
clc

Main_path = 'D:\Jazmin\MultichannelDataTanks\HIP\';
tanks = {'24_023','24_024','24_030'};
Sorter = 'KS4';

%%
for ta = 1:length(tanks)

    % Get electrode positions
    Positions = allfolders([Main_path  tanks{ta} ]);
 
    for pos = 1:length(Positions)
        clear PhyRez Spikes_blockSeparated
        tank_path = [Main_path  tanks{ta} '\' Positions{pos}];

        % Load file with details from previous combining of block data
        if exist([tank_path '\StreamSplitInfo_All.mat'])
            load([tank_path '\StreamSplitInfo_All.mat'])

            if exist([tank_path '\Sorting\' Sorter '\sorter_output\spike_times.npy'])

                % Load Phy sorted data
                clear PhyRez
                PhyRez(:,1) = readNPY([tank_path '\Sorting\' Sorter '\sorter_output\spike_times.npy']);
                PhyRez(:,2) = readNPY([tank_path '\Sorting\' Sorter '\sorter_output\spike_clusters.npy']);
                PhyRez(:,2) = PhyRez(:,2)+1; % Remove Phy zero-indexing

                % Reading Kilosort Labels
                KSLabels = readtable([tank_path '\Sorting\' Sorter '\sorter_output\cluster_KSLabel.tsv'], "FileType","text",'Delimiter', '\t');
                KSLabels.cluster_id = KSLabels.cluster_id+1;% Remove Phy zero-indexing
                % goodnums = find(contains(KSLabels.KSLabel,'good'));
                muanums = KSLabels.cluster_id(contains(KSLabels.KSLabel,'mua'));

                if ~isempty(muanums)
                    %%
                    SpikePositions = readNPY([tank_path '\Sorting\' Sorter '\sorter_output\spike_positions.npy']); % Spike shapes
                    clear neuronloc 

                    newclusters = unique(PhyRez(:,2));
                    newclusters = muanums;

                    for i = 1:length(newclusters) % For each cluster
                        spikelocs = find(PhyRez(:,2) == newclusters(i));
                        neuronloc(i) = mean(SpikePositions(spikelocs,2),'omitnan');
                     end

                    ClusterMUA = [double(newclusters)'; neuronloc];
                    save([tank_path '\ClusterGoodLocsMUA.mat'],'ClusterMUA','-v7.3')

                    disp([tanks{ta} '-' Positions{pos} '-' num2str(length(newclusters)) ' units'])


                    %% Separate spikes by which block they were recorded in

                    Ends = cumsum(StreamSplitInfo.LengthSamps); % Get end of block windows
                    Start = 1; % Setting initial start sample
                    for i = 1:length(StreamSplitInfo.Blocks)

                        curBlockName = StreamSplitInfo.Blocks{i};

                        blockWindow = [Start Ends(i)]; % Get first and last samples of this block
                        Start = Ends(i) + 1; % Setting start sample for next block

                        idx = find(PhyRez(:,1) >= blockWindow(1) & PhyRez(:,1) <= blockWindow(2));
                        Spikes_blockSeparated{i} = PhyRez(idx,:);

                    end

                    %% For each block, convert to usable data (BST format)
                    for i = 1:length(StreamSplitInfo.Blocks) % loop the different blocks

                        % Load BST
                        [bst,~,~] = bbst3(tank_path,StreamSplitInfo.Blocks{i},0); % Load bst (only epochs
                        units = muanums; % get list of unit numbers from sorted data

                        bst.Spikes = table;
                        SpikesStacked = [];
                        for ii = 1:length(units)
                            % Reduce spikes to just this unit
                            idx = find(Spikes_blockSeparated{i}(:,2) == units(ii));
                            tempSpikes = Spikes_blockSeparated{i}(idx,:);
                            tempSpikes = double(tempSpikes);

                            SpikesStacked = [SpikesStacked; tempSpikes]; %stackem
                        end

                        if i>1 %
                            SpikesStacked(:,1) = SpikesStacked(:,1) - StreamSplitInfo.LengthSamps(i-1);
                        end

                        bst.Spikes.Sample = SpikesStacked(:,1);
                        bst.Spikes.TS = SpikesStacked(:,1) / 24414.0625; % fs hard coded here
                        bst.Spikes.unit = SpikesStacked(:,2); % fs hard coded here

                        % Calculate Trial Index & raster
                        swepoff = bst.Epocs.TSOff.bind;
                        [~,~,bins] = histcounts(bst.Spikes.TS,[bst.Epocs.TSOn.bind; swepoff(end)]);
                        bst.Spikes.TrialIdx = bins;
                        bst.Spikes.TrialIdx(bst.Spikes.TrialIdx==0) = 1;
                        bst.Spikes.RasterSW = bst.Spikes.TS - bst.Epocs.TSOn.bind(bst.Spikes.TrialIdx);

                        % bst.SpikeShapes = SpikeShapesNew;

                        save([tank_path '\' StreamSplitInfo.Blocks{i} '\bst_' Sorter '_MUA.mat'],'bst')
                        disp(['Block: ' StreamSplitInfo.Blocks{i} '. Units:' num2str(length(units)) '. Spikes:' num2str(height(bst.Spikes)) '. Saved.'])
                    end
                end
            end
        end
    end
end

%%
function folders = allfolders(directory)

folders = dir(directory);
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
folders = folders(dirFlags);
folders = {folders.name};

end
