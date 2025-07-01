% Preliminary script for PCA & k-means spike sorking on single channel TDT data
% Combines multiple blocks recorded on the same neuron
% Then converts to bst format, and some example plotting
%
%  Requires TDTbin2mat and matlab SDK in matlab path, and beautify
%  function.
%
% AH 12/2023
%Modified by Jazmin Sanchez 06/2024
% All variables to change in here
close all; clear all; clc;

DATAPATH = 'D:\PRJ_HIPPOCAMPUS_JSS\Tanks\22_065';
SAVEPATH = 'D:\PRJ_HIPPOCAMPUS_JSS\Tanks\22_065';
FIGPATH = 'D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\Cortex\Figure Sorting\22_065';
SPIKESPATH = 'D:\PRJ_HIPPOCAMPUS_JSS\PCASorting\Cortex\Spike\22_065';

% Get electrode positions
files = dir(DATAPATH);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..

Blocks = {'1-2-MSC-75ms-1000ms','1-2-CASASC-75ms-1000ms','1-2-CASCDESC-75ms-1000ms','1-2-ODD-75ms-1000ms-5-6','1-2-ODD-75ms-1000ms-6-5'};

% Blocks = {'1-2-CASASC-75ms-1000ms','1-2-CASCDESC-75ms-1000ms','1-2-ODD-75ms-1000ms-8-9','1-2-ODD-75ms-1000ms-9-8'};

ThreshSTDs = 3.5; % How many standard deviations form the mean to set the threshold at.

filterFreqs = [1000 10000]; % Filter frequencies for the stream, standard is [300 5000], but [1000 to 5000] seems to remove some of the weird hippocampal fluctuations

streamdouble = [];
StreamSplitInfo = struct;
StreamSplitInfo.Blocks = Blocks;

for i = 1:length(Blocks) % For all blocks that need appending to the original binary file

    % Load TDT Block data and filter stream
    BLOCKPATH = fullfile(DATAPATH,Blocks{i});
    data = TDTbin2mat(BLOCKPATH); % read the first second of data to get the channel coun

    stream = data.streams.ECsg.data;
    streamdouble = [streamdouble double(stream)];

    StreamSplitInfo.LengthSamps(i) = length(stream);

end
fs = data.streams.ECsg.fs;
[b] = fir1(128,filterFreqs/(fs/2));
streamfilt = filtfilt(b,1,streamdouble);

% Set threshold based on stds from mean
thresh = mean(streamfilt) - (ThreshSTDs * std(streamfilt));

% Find times when the stream crosses the set threshold
idxl = streamfilt<=thresh;
idxl(1) = 0;
idx = find(idxl);
yest = streamfilt(idx-1)>thresh;
crosstimes = idx(yest);  % Final output

% remove crossings that are <1ms after a previous crossing, avoiding
% capturing a single spike twice
for iii = 2:length(crosstimes)
    if crosstimes(iii)-crosstimes(iii-1) < (0.001*fs)
        crosstimes(iii) = nan;
        %             crosstimes(iii) = [];
    end
end
crosstimes(isnan(crosstimes)) = [];

% Convert filteres spikes to matrix data with rows as individual spike waveforms
clear newdata
for ii = 1:length(crosstimes)
    if crosstimes(ii)+97 < length(streamfilt) && crosstimes(ii)-48 > 0
        newdata(ii,:) = streamfilt(crosstimes(ii)-48:crosstimes(ii)+97); % show spikes -1ms to +2ms from crossing
    else
        newdata(ii,:) = nan;
    end
end

spikewavs = newdata;

% Run PCA
[coefs,score] = pca(spikewavs);
% combine and calculate clusters of pca scores (using kmeans)
arr = score(:,1:3)';
if length(arr) > 100000
    test_arr = arr(:,1:100000)';
    msize = length(arr);
    test_arr = arr(:,randperm(msize, 100000));
else
    test_arr = arr;
end

% Kmeans clustering (unsupervised)
tiledlayout(1,4)
nexttile
clear idx C sumdist3 silh h silmean silmed
for i = 1:5 % repeat for different number of clusters
    [idx,C,sumdist3] = kmeans(test_arr',i,'Distance','sqeuclidean','Display','off','MaxIter',30000,'Replicates',100); % ,'MaxIter',1,'Replicates',1
    [silh,h] = silhouette(test_arr',idx,'sqEuclidean');
    silmean(i) = nanmean(silh);
    silmed(i) = nanmedian(silh);
end

if max(silmean) > 0.75 % if getting good separation in multiple clusters (set this manually based on testing)
    z = find(silmean == max(silmean));
    %         z = 3; % manually set number of clusters
    [idx,C,sumdist3] = kmeans(test_arr',z,'Distance','sqeuclidean','Display','off','MaxIter',30000,'Replicates',100); %
    [~,idx] = pdist2(C,arr','squaredeuclidean','Smallest',1);
    idx = idx';
else
    clear idx
    idx(1:length(arr),:) = 1;
end

% eva = evalclusters(test_arr','kmeans','CalinskiHarabasz','KList',1:10); % ,'MaxIter',1,'Replicates',1
% z = eva.OptimalK;
% [idx,C,sumdist3] = kmeans(test_arr',z,'Distance','sqeuclidean','Display','off','MaxIter',20000,'Replicates',100); %
% [~,idx] = pdist2(C,arr','squaredeuclidean','Smallest',1);
% idx = idx';

nexttile
plot(streamfilt(1:100000)) % shows first ~2s
hold on
yline(thresh)

% remove units with only a few spikes
for clu = 1:length(unique(idx))
    locs = find(idx==clu);
    if length(locs) < 10
        idx(locs) = [];
    end
end

%test the clustering on pca plot
cols = ({'r','b','g','k','m','c','y'});
clusters = unique(idx);
for i = 1:length(clusters)
    clu = clusters(i);
    locs = find(idx==clu);
    msize = numel(locs);
    nexttile(3)
    %                 biplot(coefs(:,1:3),'color','k')
    %                 hold on
    %                 biplot('Scores',score(locs(randperm(msize, 300)),1:3),'color',cols{clu})
    try
        h = biplot(coefs(:,1:3),'Scores',score(locs(randperm(msize, 300)),1:3),'color',cols{clu});
    catch
        try
            h = biplot(coefs(:,1:3),'Scores',score(locs(randperm(msize, 25)),1:3),'color',cols{clu});
        catch
            h = biplot(coefs(:,1:3),'Scores',score(locs(randperm(msize, 3)),1:3),'color',cols{clu});
        end
    end

    hold on
    %                 for k = 1:20
    %                     h(k).Color = 'k'; % Specify red as the line color
    %                     h(k).LineStyle = '--'; % Specify red as the line color
    %                 end
    clulocs{clu} = find(idx==clu);
    nexttile(4)
    plot(mean(spikewavs(clulocs{clu},:)),cols{clu},'DisplayName',num2str(length(clulocs{clu})))
    hold on

end
legend
box off
beautify
cd (FIGPATH)
saveas(gcf,'Sorting.fig')%save fig

%% Jazmin sorting part 2, convert to bst format

for i = 1:length(clusters)

    %load in spike times and clusters
    % find locations for cluster
    locs = find(idx == clusters(i));

    clear Spikes_all
    Spikes_all(:,1) = crosstimes(locs);% load the spike times
    Spikes_all(:,2) = clusters(i);% load the spike cluster

    % Create bst for this block
    save_dir = DATAPATH;

    % Create superblocks
    superblocks = build_sb_jazmin(DATAPATH, save_dir, StreamSplitInfo.Blocks(1),'dontsave');
    % Load BST
    base_bst = bbst_jazmin(superblocks,save_dir,StreamSplitInfo.Blocks(1),1); % Load data for channel 1 - just getting epoks not spikes so could use any channel
    base_bst = rmfield(base_bst, 'Spikes'); % Remove TDT-created spikes

    Ends = cumsum(StreamSplitInfo.LengthSamps); % Get end of block windows

    if length(StreamSplitInfo.Blocks) > 1
        % For each block, convert to usable data (BST format)
        for ii = 2:length(StreamSplitInfo.Blocks) % loop the different blocks

            % Create superblocks
            superblocks_temp = build_sb_jazmin(DATAPATH, save_dir, StreamSplitInfo.Blocks(ii),'dontsave');
            % Load BST
            append_bst = bbst_jazmin(superblocks_temp,save_dir,StreamSplitInfo.Blocks(ii),1); % Load data for channel 1 - just getting epoks not spikes so could use any channel
            append_bst = rmfield(append_bst, 'Spikes'); % Remove TDT-created spikes

            % append epochs.values
            base_bst.Epocs.Values = [base_bst.Epocs.Values ; append_bst.Epocs.Values];

            append_bst.Epocs.TSOn.levl = append_bst.Epocs.TSOn.levl + (Ends(ii)/fs);
            append_bst.Epocs.TSOff.levl = append_bst.Epocs.TSOff.levl + (Ends(ii)/fs);

            append_bst.Epocs.TSOn.freq = append_bst.Epocs.TSOn.freq + (Ends(ii)/fs);
            append_bst.Epocs.TSOff.freq = append_bst.Epocs.TSOff.freq + (Ends(ii)/fs);

            base_bst.Epocs.TSOn = [base_bst.Epocs.TSOn ; append_bst.Epocs.TSOn];
            base_bst.Epocs.TSOff = [base_bst.Epocs.TSOff ; append_bst.Epocs.Values];

        end
    end

    bst = base_bst; %

    % Reduce spikes to just this unit
    tempSpikes = Spikes_all(:,1);

    bst.Spikes = table;
    bst.Spikes.Sample = tempSpikes(:,1);
    bst.Spikes.TS = tempSpikes(:,1) / fs;

    % Calculate TrialIdx & Raster timings
    for sp = 1:length(bst.Spikes.TS) % for each spike
        bst.Spikes.TrialIdx(sp) = length(find(bst.Spikes.TS(sp) > bst.Epocs.TSOn.freq));
        if bst.Spikes.TrialIdx(sp) > 0
            bst.Spikes.RasterSW(sp) = bst.Spikes.TS(sp) - bst.Epocs.TSOn.freq(bst.Spikes.TrialIdx(sp));
        end
    end

    bst = BST_AddType_jazmin(bst);

    save([save_dir '\bst_cluster_' num2str(i)],'bst','-v7.3')

    % Select trials and save data table
    Epoks = bst.Epocs.Values;

    % Get the F1 and F2 used
    SelectedTrials = BST_TS(bst,'type','ODD-DES-DEV');
    F1 = unique(Epoks.freq(SelectedTrials(1)));
    SelectedTrials = BST_TS(bst,'type','ODD-ASC-DEV');
    F2 = unique(Epoks.freq(SelectedTrials(1)));

    sweeplength = bst.Epocs.TSOn.freq(21) - bst.Epocs.TSOn.freq(20); % Get the sweep length in seconds

    % Create a struct of spike data for each condition type
    spikedata = table;
    conditions = uniqueCellVector(Epoks.type);
    conditions(cellfun(@isempty,conditions))=[];

    for co = 1:length(conditions)
        SelectedTrials = BST_TS(bst,'type',conditions{co});
        SpikeTimes = 1000 * BST_GS(bst,SelectedTrials,'SW');
        spikedata.(conditions{co}) = {SpikeTimes};
    end

    % Save the spikedata structure
    cd (SPIKESPATH)
    save(['spikedata' num2str(i) '.mat'],"spikedata")

    % Plot histograms
    devspikes = [spikedata.('ODD-ASC-DEV'){1}; spikedata.('ODD-DES-DEV'){1}];
    stdspikes = [spikedata.('ODD-ASC-STD'){1}; spikedata.('ODD-DES-STD'){1}];
    ctrspikes = [spikedata.('CASC-ASC-F1'){1}; spikedata.('CASC-ASC-F2'){1}; spikedata.('CASC-DES-F1'){1}; spikedata.('CASC-DES-F2'){1}];

    histwidth = 25;
    histedges = 0:histwidth:1000;
    histcenters = histedges(1:end-1)+histwidth/2;

    tiledlayout(1,3)
    nexttile
    histdat = histcounts(stdspikes,histedges);
    histdat = (histdat/80)/(histwidth/40); % normalise to Hz
    bar(histcenters,histdat,'k','BarWidth',1)
    hold on
    xticks(histedges(1:5:end))
    title('STD')
    temp = get(gca,'ylim');
    ydat(1) = temp(2);

    nexttile
    histdat = histcounts(ctrspikes,histedges);
    histdat = (histdat/160)/(histwidth/40); % normalise to Hz
    bar(histcenters,histdat,'k','BarWidth',1)
    hold on
    xticks(histedges(1:5:end))
    title('CTR')
    temp = get(gca,'ylim');
    ydat(2) = temp(2);

    nexttile
    histdat = histcounts(devspikes,histedges);
    histdat = (histdat/80)/(histwidth/40); % normalise to Hz
    bar(histcenters,histdat,'k','BarWidth',1)
    hold on
    xticks(histedges(1:5:end))
    title('DEV')
    temp = get(gca,'ylim');
    ydat(3) = temp(2);

    nexttile(1)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(2)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(3)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')

    beautify
    cd (FIGPATH)
    saveas(gcf,['Histograms' num2str(i) '.pdf'])
    %save fig
end


