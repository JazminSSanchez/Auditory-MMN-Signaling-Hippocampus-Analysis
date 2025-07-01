clear all 
load('D:\Jazmin\MultichannelDataTanks\HIP\24_023\1\1-1-STIM\bst_KS4.mat')
% Created by Jazmin Sanchez 2025
% Step 1: Extract spike info for unit 83
n = 1;
for i = 1:size(bst.Spikes.RasterSW,1)
    if bst.Spikes.unit(i) == 83
        Spikes_Tr_ms(n,1) = bst.Spikes.TrialIdx(i); % Trial index
        Spikes_Tr_ms(n,2) = bst.Spikes.RasterSW(i); % spike time
        Spikes_Tr_ms(n,3) = bst.Spikes.TS(i);       % spike timestamp
        n = n + 1;
    end
end

% Step 2: Identify trial types
ODDRTrials = 461:420;

STDtrials = [];
DEVtrials = [];
m = 1;
k = 1;
for j = ODDRTrials
    if contains(bst.Epocs.Values.type{j}, 'STD')
        STDtrials(m,1) = j;
        m = m + 1;
    elseif contains(bst.Epocs.Values.type{j}, 'DEV')
        DEVtrials(k,1) = j;
        k = k + 1;
    end
end

% Step 3: Create raster plot
figure; hold on;

% Reverse y-axis mapping: trial 400 -> y=1, trial 0 -> y=401
reverseTrialMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
for t = ODDRTrials
    reverseTrialMap(t) = 420 - t + 1;
end

% Step 4: Plot spikes as dots
for i = 1:size(Spikes_Tr_ms,1)
    trialNum = Spikes_Tr_ms(i,1);
    spikeTime = Spikes_Tr_ms(i,2);

    if ismember(trialNum, ODDRTrials)
        y = reverseTrialMap(trialNum);
        if ismember(trialNum, STDtrials)
            plot(spikeTime, y, 'b.', 'MarkerSize', 10); % STD: blue
        elseif ismember(trialNum, DEVtrials)
            plot(spikeTime, y, 'r.', 'MarkerSize', 10); % DEV: red
        else
            plot(spikeTime, y, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 10); % Other: gray-ish
        end
    end
end

% Step 5: Formatting
xlabel('Time (s)');
ylabel('Trial (400 to 0)');
title('Spike Raster Plot: STD (blue), DEV (red), Other (gray)');
set(gca, 'YDir', 'normal'); % Important: reverse Y axis (trial 400 on top)
ylim([1, length(ODDRTrials)]); % Keep full trial range visible


cd 'D:\Jazmin\MultichannelDataTanks\HIP\DotRaster';

exportgraphics(gcf, 'raster_plot_83.pdf', 'ContentType', 'vector');
