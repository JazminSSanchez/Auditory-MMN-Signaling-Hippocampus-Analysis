%% LFP and CSD data saving
% Author: Adam Hockley, June 2023
% Changed, commented and cleaned by Jazmin Sanchez Sanchez, January 2025
% Description:
%   Load TDT tank data for local field potential (LFP) and current source density (CSD) analysis.
%   Saves the processed data as .mat files for downstream plotting tools.

%% Configuration
tanks = {'24_023','24_024','24_030'};  % List of tank IDs
saveCSDdat = true;     % Flag to save CSD data
saveLFPdat = true;     % Flag to save LFP data

%% Process each tank
for ta = 1:numel(tanks)
    tankPath = fullfile('D:\Jazmin\MultichannelDataTanks\HIP\', tanks{ta}, filesep);
    savePath = ('D:\Jazmin\MultichannelDataTanks\HIP\testJune2025')
    electrodeDirs = listSubfolders(tankPath);

    for el = 1:numel(electrodeDirs)
        blocks = listSubfolders(fullfile(tankPath, electrodeDirs{el}));
        % Only process STIM blocks
        blocks = blocks(contains(blocks, 'STIM'));

        for bl = 1:numel(blocks)
            blockPath = fullfile(tankPath, electrodeDirs{el}, blocks{bl});
            % Load streaming and epoc data (streams type 2 & 4)
            data = TDTbin2mat(blockPath, 'TYPE', [2 4]);
            
            % Build bst structure with types assigned
            bst = bbst3(fullfile(tankPath, electrodeDirs{el}), blocks{bl}, 0);
            bst = BST_AddType(bst);

            % Extract event onset indices by condition
            locs = extractEpocLocations(bst);

            %% Compute and save CSD data
            if saveCSDdat
                CSDdat = computeCSD(data, locs.dev, 39);
                safeSave(blockPath, 'csd_data', CSDdat);
            end

            %% Compute and save LFP data per condition
            if saveLFPdat
                lfpTable = computeLFP(data, locs, 39);
                safeSave(blockPath, 'lfp_data', lfpTable);
            end
        end
    end
end

%% Utility functions

function folders = listSubfolders(directory)
    % Return list of valid subdirectory names
    d = dir(directory);
    flags = [d.isdir] & ~ismember({d.name}, {'.','..'});
    folders = {d(flags).name};
end

function locs = extractEpocLocations(bst)
    % Create struct mapping event types to epoc indices
    types = bst.Epocs.Values.type;
    locs.dev = find(strcmp(types, 'ODD-ASC-DEV') | strcmp(types, 'ODD-DES-DEV'));
    locs.devPdev = find(strcmp(types, 'ODD-ASC-DEV-P') | strcmp(types, 'ODD-DES-DEV-P'));
    locs.ctr = find(strcmp(types, 'MS-F1') | strcmp(types, 'MS-F2'));
    locs.cascasc = find(strcmp(types, 'CASC-ASC-F1') | strcmp(types, 'CASC-ASC-F2'));
    locs.cascdes = find(strcmp(types, 'CASC-DES-F1') | strcmp(types, 'CASC-DES-F2'));

    % Handle omission F1/F2 events
    omF1 = find(strcmp(types, 'OM-F1'));
    omF2 = find(strcmp(types, 'OM-F2'));
    if ~isempty(omF1)
        locs.OmF1Ran = omF1(1:40);
        locs.OmF1Per = [omF1(41:end), omF1(end)+(omF1(end)-omF1(end-1))];
        locs.OmF2Ran = omF2(1:40);
        locs.OmF2Per = [omF2(41:end), omF2(end)+(omF2(end)-omF2(end-1))];
    end
end

function CSDdat = computeCSD(data, devIndices, downsample)
    % Filter onset-related LFP data and compute average CSD
    fs = data.streams.SU_2.fs;
    [b, a] = butter(2, [0.1 30] / (fs/2), 'bandpass');
    onsets = round(data.epocs.Us1_.onset(devIndices) * fs);
    numCh = size(data.streams.SU_2.data, 1);
    nOnsets = numel(onsets);
    snippetLen = round(fs) * 0.5;  % half second before and after

    sweepdat = nan(nOnsets, 2*snippetLen/downsample);
    CSDdat = nan(2*snippetLen/downsample, numCh);

    for ch = 1:numCh
        trace = double(data.streams.SU_2.data(ch, 1:onsets(end)+snippetLen));
        sweeps = arrayfun(@(o) decimate(filtfilt(b, a, trace(o-snippetLen:o+snippetLen)), downsample), ...
                          onsets, 'UniformOutput', false);
        sweeps = cat(1, sweeps{:});
        CSDdat(:, ch) = mean(sweeps, 1)';
        fprintf('Channel %d processed for CSD\n', ch);
    end
end

function lfpTable = computeLFP(data, locs, downsample)
    % Compute per-condition LFP snippets for each channel
    fs = data.streams.SU_2.fs;
    [b, a] = butter(2, [0.1 30] / (fs/2), 'bandpass');
    conditions = fieldnames(locs);
    numCh = size(data.streams.SU_2.data, 1);

    lfpTable = table;
    for ch = 1:numCh
        trace = double(data.streams.SU_2.data(ch, :));
        for c = 1:numel(conditions)
            idxs = locs.(conditions{c});
            onsets = round(data.epocs.Us1_.onset(idxs) * fs);
            snippets = arrayfun(@(o) decimate(filtfilt(b, a, trace(max(1,o-fs):min(end,o+2*fs))), downsample), ...
                                onsets, 'UniformOutput', false);
            lfpTable.(conditions{c}){ch} = cat(1, snippets{:});
        end
        fprintf('Channel %d processed for LFP\n', ch);
    end
end

function safeSave(path, filename, dataVar)
    % Saves dataVar into file under path\filename.mat
    if ~exist(path, 'dir')
        mkdir(path);
    end
    save(fullfile(path, [filename '.mat']), inputname(3), '-v7.3');
end