function bst = bbst_jazmin(superblocks,tank_path,block_list,ch)

% AH2019 - largely undocumented major improvement to superspiketrain.
bst = struct;
bst.tank = tank_path;
bst.Block = block_list;
bst.channel = ch;

bst.EpocNames = {};
bst.Epocs = struct;

% Get just the spikes for this channel
toDelete = superblocks{1}.chan ~= ch;
superblocks{1}(toDelete,:) = [];
bst.Spikes = superblocks{1};

for bl = 1:length(block_list) % loop blocks

    epochs = TDTbin2mat([tank_path '/' block_list{bl}],'Type',[2],'Verbose',0); % Load epochs
    fns = fieldnames(epochs.epocs);
    nsweps = length(epochs.epocs.Freq.data);
    if bl == 1 % if it's the first block, add epoc names
        bst.EpocNames{1} = 'bind';
        for ie = 1:length(fns)
            bst.EpocNames{end+1,1} = lower(fns{ie});
        end
        bst.Epocs.Values = table;
        bst.Epocs.TSOn = table;
        bst.Epocs.TSOff = table;
    end
    
    h = height(bst.Epocs.Values);
    clear blocks bind
%     blocks(1:nsweps) = block_list{bl};
    bind(1:nsweps) = bl;
    
    warning off
    % Add data values to table
    bst.Epocs.Values.Block(h+1:h+nsweps) = block_list(bl);
    bst.Epocs.Values.bind(h+1:h+nsweps) = bind';
    for fi = 1:length(fns)
        tempdata = epochs.epocs.(fns{fi}).data;
        bst.Epocs.Values.(lower(fns{fi}))(h+1:h+length(tempdata)) = tempdata;

    end
    
    % Add TSOn values to table
    bst.Epocs.TSOn.Block(h+1:h+nsweps) = block_list(bl);
    bst.Epocs.TSOn.bind(h+1:h+nsweps) = epochs.epocs.Freq.onset(1:nsweps); % 2 to end?
    for fi = 1:length(fns)
        tempdata = epochs.epocs.(fns{fi}).onset;
        bst.Epocs.TSOn.(lower(fns{fi}))(h+1:h+length(tempdata)) = tempdata;

    end

    % Add TSoff values to table
    bst.Epocs.TSOff.Block(h+1:h+nsweps) = block_list(bl);
    bst.Epocs.TSOff.bind(h+1:h+nsweps) = epochs.epocs.Freq.offset(1:nsweps); % 2 to end?
    for fi = 1:length(fns)
        tempdata = epochs.epocs.(fns{fi}).offset;
        bst.Epocs.TSOff.(lower(fns{fi}))(h+1:h+length(tempdata)) = tempdata;

    end
    warning on
    
end

bst.NTrials = height(bst.Epocs.Values);
bst.Spikes.Properties.VariableNames(4) = {'TS'};
bst.Spikes.Properties.VariableNames(5) = {'SortCodes'};

%% Add raster data to bst

%%% Assign timestamps to trials
% Loop through each timestamp, find which trial it belongs to
% and subtract that swep timestamp
clear TrialCode
TrialCode = zeros(length(bst.Spikes.TS),1);

for bl = 1:length(bst.Block)
    
    %Determine off time in each trial for current block
    swepoff = bst.Epocs.TSOff.freq(ismember(bst.Epocs.TSOff.Block,bst.Block(bl)));
    
    %Uses histogram counting function to seperate timestamp
    %data into bins
    [~,~,bins] = histcounts(bst.Spikes.TS(bst.Spikes.BlockIdx==bl),[bst.Epocs.TSOn.freq(ismember(bst.Epocs.TSOn.Block,bst.Block(bl))); swepoff(end)]);
    
    %assign here the trial number that each spike belongs to
    TrialCode(bst.Spikes.BlockIdx==bl) = bins+(bins~=0).*(find(ismember(bst.Epocs.Values.Block,bst.Block(bl)),1)-1);
    
end

bst.Spikes.TrialIdx = TrialCode;
bst.Spikes(TrialCode==0,:) = [];

%%% Convert timestamps to rasters
bst.Spikes.RasterSW = bst.Spikes.TS - bst.Epocs.TSOn.freq(bst.Spikes.TrialIdx);
if ismember('s1ig',fieldnames(bst.Epocs.TSOn))
    bst.Spikes.RasterS1=bst.Spikes.TS-bst.Epocs.TSOn.s1ig(bst.Spikes.TrialIdx);
elseif ismember('s1on',fieldnames(bst.Epocs.TSOn))
    bst.Spikes.RasterS1=bst.Spikes.TS-bst.Epocs.TSOn.s1on(bst.Spikes.TrialIdx);
end
if ismember('s2ig',fieldnames(bst.Epocs.TSOn))
    bst.Spikes.RasterS2=bst.Spikes.TS-bst.Epocs.TSOn.s2ig(bst.Spikes.TrialIdx);
elseif ismember('s2on',fieldnames(bst.Epocs.TSOn))
    bst.Spikes.RasterS2=bst.Spikes.TS-bst.Epocs.TSOn.s2on(bst.Spikes.TrialIdx);
end
if ismember('etyp',fieldnames(bst.Epocs.TSOn))
    bst.Spikes.RasterE=bst.Spikes.TS-bst.Epocs.TSOn.etyp(bst.Spikes.TrialIdx);
end


nspikesloaded = size(bst.Spikes,1);
if nspikesloaded < 100
    warning('N spikes is < 100');
elseif nspikesloaded > 1000000
    warning('N spikes is > 1000000');
end


end