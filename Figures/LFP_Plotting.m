%% Separate LFP Plotting Preparation
% Author: Adam Hockley, June 2023
% Cleaned, commented, and adapted by Jazmin Sanchez Sanchez, January 2025
%
% Description:
%   This script loads LFP data from multiple experiments,
%   applies common median referencing, and identifies noisy channels
%   based on standard deviation thresholds.
%   To look at plots or stats, it is necesarry to add debuggin stops 


%% Setup
clear;

% List of tank identifiers to process
tanks = {'24_023', '24_024', '24_030'};

% Layer options (not used in current script but could be used for layer-specific plotting)
layer = {'All'};

% Optional threshold for data cutoff or artifact rejection (currently unused)
cutoff = -0.00007;

%% Loop through each tank
for ta = 1:length(tanks)

    tankPath = fullfile('D:\Jazmin\MultichannelDataTanks\HIP', tanks{ta});
    folders = dir(tankPath);
    
    % Get electrode folder names (skip '.' and '..')
    electrodeLocs = {folders(3:end).name};

    %% Loop through each electrode folder
    for el = 1:length(electrodeLocs)

        electrodePath = fullfile(tankPath, electrodeLocs{el});
        folders = dir(electrodePath);
        blocks = {folders(3:end).name};

        % Keep only STIM blocks
        blocks = blocks(contains(blocks, 'STIM'));

        %% Loop through each STIM block
        for bl = 1:length(blocks)

            blockPath = fullfile(electrodePath, blocks{bl});

            % Load precomputed CSD data
            load(fullfile(blockPath, 'csd_data.mat'), 'CSDdat');

            %% Common Median Referencing
            % Subtract the median across channels at each time point
            for ch = 1:width(CSDdat)
                CSDdat(:, ch) = CSDdat(:, ch) - nanmedian(CSDdat, 2);
            end

            %% Detect Noisy Channels
            % Compute standard deviation across a baseline window (first 100 samples)
            stddevs = std(CSDdat(1:100, :), 0, 1);

            % Identify channels with std dev greater than the 95th percentile
            dchans = find(stddevs > prctile(stddevs, 95));

        end
    end
end

%%
for la = 1:length(layer)

    count = 0;
    for ta = 1:length(tanks)

        tankpath = ['D:\Jazmin\MultichannelDataTanks\HIP\' tanks{ta} '\'];
        electrodelocs = allfolders(tankpath);

        if exist('blocksAllManual')
            electrodelocs = blocksAllManual{ta};  % or manual block selection
        else
            % Auto All
            folders = dir(tankpath);
            electrodelocs = {folders(3:end).name};
        end

        for el = 1:length(electrodelocs)

            blocks = allfolders([tankpath electrodelocs{el}]);
            blocks = blocks(contains(blocks,'STIM'));
            blocks = blocks(~contains(blocks,'1000'));

            for bl = 1:length(blocks)

                dirpath = [tankpath electrodelocs{el} '\' blocks{bl}];
                load([dirpath '\lfp_data'])
               channels = [1:32];
               channels = channels(find(channels~=dchans(1)));
               if size(dchans,2) == 2
                   channels = channels(find(channels~=dchans(2)));
               end

               if ta==3 && el==1
                   channels (5) = [];
                   channels (24) = [];
               end

                if ~isempty(channels)
                    clear tempdev tempdevP tempdevPstd tempctr tempctrled
                    for ch = 1:length(channels)

                        channel = channels(ch);

                        %% Get the normal way
                        % TESTING WITH MEDIAN
                        tempdev(ch,:) = mean(lfpdat.dev{channel});
                        tempdevP(ch,:) = mean(lfpdat.devPdev{channel});
                        tempctr(ch,:) = mean(lfpdat.ctr{channel});
                        tempdevP(ch,:) = mean(lfpdat.devPdev{channel});


                    end
                    %% Graph of channels with neurons with significant CSI
                    load('D:\Jazmin\MultichannelDataTanks\HIP\SignificantChannels.mat')
                    N_24_023_2 = [4,7,2,25,26,28,29];
                    N_24_023_3 = [12,30];
                    N_24_024_3 = [2];
                    N_24_024_5 = [29];
                    N_24_024_6 = [14];
                    N_24_030_6 = [3];


                    if ta==1 && el==1
                        devM_S_DG = tempdev(N_24_023_2,:);
                        devMP_S_DG = tempdevP(N_24_023_2,:);
                        ctrM_S_DG = tempctr(N_24_023_2,:);
                    elseif ta==1 && el==2
                        devM_S_DG = [devM_S_DG;tempdev(N_24_023_3,:)];
                        devMP_S_DG = [devMP_S_DG;tempdevP(N_24_023_3,:)];
                        ctrM_S_DG = [ctrM_S_DG;tempctr(N_24_023_3,:)];
                    elseif ta==2 && el==1
                        devM_S_DG = [devM_S_DG;tempdev(N_24_024_3,:)];
                        devMP_S_DG = [devMP_S_DG;tempdevP(N_24_024_3,:)];
                        ctrM_S_DG = [ctrM_S_DG;tempctr(N_24_024_3,:)];
                    elseif ta==2 && el==2
                        devM_S_DG = [devM_S_DG;tempdev(N_24_024_5,:)];
                        devMP_S_DG = [devMP_S_DG;tempdevP(N_24_024_5,:)];
                        ctrM_S_DG = [ctrM_S_DG;tempctr(N_24_024_5,:)];
                    end
                    if ta==2 && el==3
                        devM_S_DG = [devM_S_DG;tempdev(N_24_024_6,:)];
                        devMP_S_DG = [devMP_S_DG;tempdevP(N_24_024_6,:)];
                        ctrM_S_DG = [ctrM_S_DG;tempctr(N_24_024_6,:)];
                    elseif ta==3 && el==3
                        devM_S_DG = [devM_S_DG;tempdev(N_24_030_6,:)];
                        devMP_S_DG = [devMP_S_DG;tempdevP(N_24_030_6,:)];
                        ctrM_S_DG = [ctrM_S_DG;tempctr(N_24_030_6,:)];
                    end
                    %% 
%% Graph of channels with neurons with significant CSI
                    load('D:\Jazmin\MultichannelDataTanks\HIP\SignificantChannels.mat')
                    N_24_023_2 = [10,14,16];

                    if ta==1 && el==1
                        devM_S_CA1 = tempdev(N_24_023_2,:);
                        devMP_S_CA1 = tempdevP(N_24_023_2,:);
                        ctrM_S_CA1 = tempctr(N_24_023_2,:);
                    end
                        %% 
                    % Remove channels with very small LFPs (in noise floor)
                    RemAmps = peak2peak(tempdev(783:1033)');
                    RemIdx = RemAmps < cutoff; % approx 0.0002 is good

                    tempdev(RemIdx,:) = [];
                    tempdevP(RemIdx,:) = [];
                    tempctr(RemIdx,:) = [];

                    count = count+1;

                    % Get av waveforms
                   
                    if ta==1 && el==1
                        devM = tempdev;
                        devMP = tempdevP;
                        ctrM = tempctr;
                    else
                        devM = [devM;tempdev];
                        devMP = [devMP;tempdevP];
                        ctrM = [ctrM;tempctr];
                    end

                    dev(count,:) = mean(tempdev);
                    devP(count,:) = mean(tempdevP);
                    % devledstd(count,:) = mean(tempdevPstd);
                    ctr(count,:) = mean(tempctr);
                    % ctrled(count,:) = mean(tempctrled);
                end
            end
        end
    end

                    %% Added plotting for mean of that layer for each electrode location
                    timevector = (1:length(tempdev))/0.626;
                    timevector = timevector-1000;

                    stdSamps = 1:626;
                    devSamps = 626:1252;
                    postSamps = 1252:length(timevector);
            figure('Position', get(0, 'Screensize'));
            tempdevT = tempdev';
            SEM_m = mean(mean(tempdev))/sqrt(length(channels));
            
            savpath = [dirpath '\LFP'];
            
                        %// Define your data
                        data = tempdevT';
                        
                        %// Define integer grid of coordinates for the above data
                        [X,Y] = meshgrid(1:size(data,2), 1:size(data,1));
                        
                        %// Define a finer grid of points
                        [X2,Y2] = meshgrid(1:0.01:size(data,2), 1:0.01:size(data,1));
                        
                        %// Interpolate the data and show the output
                        outData = interp2(X, Y, data, X2, Y2, 'linear');
                        imagesc(outData);
                        
                        %// Cosmetic changes for the axes
                        set(gca, 'XTick', linspace(1,size(X2,2),size(X,2))); 
                        set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
                        set(gca, 'XTickLabel', 1:size(X,2));
                        set(gca, 'YTickLabel', 1:size(X,1));
                        
                        %// Add colour bar
                        colorbar;
                        print(gcf,'-vector','-depsc',[savpath '\heatmap' '.eps'])
                        print(gcf,'-dpng','-r150',[savpath '\heatmap' '.png'])
                        print(gcf,'-dpng','-r150',[savpath '\heatmap' '.fig'])
                        close all

            figure
            StacjedLFP1_15 = stackedplot(timevector,tempdevT(:,1:15));

            for i = 1:15
                StacjedLFP1_15.AxesProperties(i).YLimits = [min(min(tempdev)),max(max(tempdev))];
            end

            if ~exist(savpath,'dir'), mkdir(savpath); 
            end

            print(gcf,'-vector','-depsc',[savpath '\StakedLFP_BP4_0.1-30Hz_1' '.eps'])
            print(gcf,'-dpng','-r150',[savpath '\StakedLFP_BP4_0.1-30Hz_1' '.png'])
            print(gcf,'-dpng','-r150',[savpath '\StakedLFP_BP4_0.1-30Hz_1' '.fig'])
            close all

            figure('Position', get(0, 'Screensize'));
            StacjedLFP16_end = stackedplot(timevector,tempdevT(:,16:length(channels)));

            for i = 1:length(channels)-15
                StacjedLFP16_end.AxesProperties(i).YLimits = [min(min(tempdev)),max(max(tempdev))];
            end

            print(gcf,'-vector','-depsc',[savpath '\StakedLFP_BP4_0.1-30Hz_2' '.eps'])
            print(gcf,'-dpng','-r150',[savpath '\StakedLFP_BP4_0.1-30Hz_2' '.png'])
            print(gcf,'-dpng','-r150',[savpath '\StakedLFP_BP4_0.1-30Hz_2' '.fig'])
            close all

            %% 
                    figure
                    tiledlayout(1,4)
                    lims = [];
                    nexttile
                    plot(timevector(stdSamps),mean(tempdev(:,stdSamps))*1000,'b')
                    hold on
                    plot(timevector(devSamps),mean(tempdev(:,devSamps))*1000,'r')
                    plot(timevector(postSamps),mean(tempdev(:,postSamps))*1000,'k')

                    plot(timevector(stdSamps),mean(tempdevP(:,stdSamps))*1000,'b--')
                    plot(timevector(devSamps),mean(tempdevP(:,devSamps))*1000,'r--')
                    plot(timevector(postSamps),mean(tempdevP(:,postSamps))*1000,'k--')

                    lims = ylim;
                    ylim([-0.06 0.08])

                    nexttile
                    plot(timevector,mean(tempctr)*1000,'g')
                    hold on

                       ylim([-0.06 0.08])
                    %% getting bar data PEAK 2 PEAK
                    ampsstd = peak2peak(tempdev(:,156:406)');
                    ampsstdP = peak2peak(tempdevP(:,156:406)');

                    ampsdev = peak2peak(tempdev(:,783:1033)');
                    ampsdevP = peak2peak(tempdevP(:,783:1033)');

                    ampspost = peak2peak(tempdev(:,1409:1659)');
                    ampspostP = peak2peak(tempdevP(:,1409:1659)');

                    ampsctr = peak2peak(tempctr(:,156:406)');

                    nexttile
                    x=categorical({'STD';'DEV';'CTR'});
                    x=reordercats(x,{'STD';'DEV';'CTR'});
                    y=[mean(ampsstd) mean(ampsstdP);  mean(ampsdev) mean(ampsdevP);  mean(ampsctr) 0] *1000;
                    errorplus=([std(ampsstd) std(ampsstdP);  std(ampsdev) std(ampsdevP);  std(ampsctr) 0] *1000 ) / sqrt(length(ampsctr));
                    errorminus=errorplus;
                    b = bar(x, y, 0.8, 'FaceColor' , 'flat');
                    clear ctrl ydt
                    for k1 = 1:size(y,2)
                        ctrl(k1,:) = bsxfun(@plus, [1 2 3], b(k1).XOffset');
                        ydt(k1,:) = b(k1).YData;
                    end
                    hold on
                    errorbar(ctrl, ydt, errorplus', '.k')

                     % 1/3
                    [p, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[2]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis([ampspost; ampspostP]' ,[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end


                    figtexts('','LFP Amplitude (mV)','')
                    b(1).CData(1,:) = [0 0.4470 0.7410];
                    b(2).CData(1,:) = [0 0.4470 0.7410]/0.7;
                    b(1).CData(2,:) = [0.8500 0.3250 0.0980];
                    b(2).CData(2,:) = [0.8500 0.3250 0.0980]/0.7;
                    b(1).CData(3,:) = [0.4660 0.6740 0.1880];
                    b(2).CData(3,:) = [0.4660 0.6740 0.1880]/0.7;

                    %% Now PE plotting (new way)

                    euclidnorm = sqrt((ampsdev.^2)+(ampsstd.^2)+(ampsctr.^2));
                    euclidnormP = sqrt((ampsdevP.^2)+(ampsstdP.^2)+(ampsctr.^2));

                    ampsdev = ampsdev ./ euclidnorm;
                    ampsctr = ampsctr ./ euclidnorm;
                    ampsstd = ampsstd ./ euclidnorm;
                    ampsdevP = ampsdevP ./ euclidnormP;
                    ampsstdP = ampsstdP ./ euclidnormP;

                    iPE = ampsdev - ampsctr;
                    iRS = ampsctr - ampsstd;
                    iPEp = ampsdevP - ampsctr;
                    iRSp = ampsctr - ampsstdP;

                    iMM = iPE + iRS;
                    iMMp = iPEp + iRSp;

                    nexttile
                    violindata = [iPE', iPEp', iRS', iRSp', iMM', iMMp'];

                    [p, tbl, stats] = kruskalwallis(violindata,[],'off');
                    data = multcompare(stats,'display','off');

                    vp = violinplot(violindata, {'iPE','iPE P','iRS','iRS P','iMM','iMM P'});
                    vpColors = [[255 153 0]*0.85; [255 153 0]; [0 153 204]*0.85; [0 153 204]; [229 0 126]*0.85; [229 0 126]];
                    vpColors = vpColors/255;
                    for bw = 1:length(vp)
                        vp(bw).BoxWidth = 0.05;
                        vp(bw).WhiskerPlot.LineWidth = 2;
                        vp(bw).ViolinColor{1} = vpColors(bw,:);
                    end

                    % 1/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,1:2),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1 2]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,3:4),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[3 4]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,5:6),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','Mismatch index','')

                    %%
                    sgtitle(['Channels = ' num2str(length(iPE))])

                    beautify

                    savpath = [dirpath '\LFP'];
                    if ~exist(savpath,'dir'), mkdir(savpath); end

                    %                     exportgraphics(gcf,[savpath '\_' layer{la} 'NEW.eps'])% x = imresize(x,0.4811161);
                    print(gcf,'-vector','-depsc',[savpath '\LFP_BP4_0.1-30Hz' '.eps'])
                    print(gcf,'-dpng','-r150',[savpath '\LFP_BP4_0.1-30Hz' '.png'])
                    close all

                    %% MAtrix analysis of significant reponses
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Matrix created before around line 173
                        % devM;
                        % devMP;
                        % ctrM;
                    %% Added plotting for mean of that layer for each electrode location
                    timevector = (1:length(devM))/0.626;
                    timevector = timevector-1000;

                    stdSamps = 1:626;
                    devSamps = 627:1252;
                    postSamps = 1253:length(timevector);

                    figure
                    tiledlayout(1,4)
                    lims = [];
                    nexttile
                    plot(timevector(stdSamps),mean(devM(:,stdSamps))*1000,'b')
                    hold on
                    plot(timevector(devSamps),mean(devM(:,devSamps))*1000,'r')
                    plot(timevector(postSamps),mean(devM(:,postSamps))*1000,'k')

                    plot(timevector(stdSamps),mean(devMP(:,stdSamps))*1000,'b--')
                    plot(timevector(devSamps),mean(devMP(:,devSamps))*1000,'r--')
                    plot(timevector(postSamps),mean(devMP(:,postSamps))*1000,'k--')

                    lims = ylim;
                     ylim([-0.06 0.08])

                    nexttile
                    plot(timevector,mean(ctrM)*1000,'g')
                    hold on
                    ylim([-0.06 0.08])

                    %% getting bar data PEAK 2 PEAK
                    ampsstd = peak2peak(devM(:,156:406)');
                    ampsstdP = peak2peak(devMP(:,156:406)');

                    ampsdev = peak2peak(devM(:,783:1033)');
                    ampsdevP = peak2peak(devMP(:,783:1033)');

                    ampspost = peak2peak(devM(:,1409:1659)');
                    ampspostP = peak2peak(devMP(:,1409:1659)');

                    ampsctr = peak2peak(ctrM(:,156:406)');

                    nexttile
                    x=categorical({'STD';'DEV';'CTR'});
                    x=reordercats(x,{'STD';'DEV';'CTR'});
                    y=[mean(ampsstd) mean(ampsstdP);  mean(ampsdev) mean(ampsdevP);  mean(ampsctr) 0] *1000;
                    errorplus=([std(ampsstd) std(ampsstdP);  std(ampsdev) std(ampsdevP);  std(ampsctr) 0] *1000 ) / sqrt(length(ampsctr));
                    errorminus=errorplus;
                    b = bar(x, y, 0.8, 'FaceColor' , 'flat');
                    clear ctrl ydt
                    for k1 = 1:size(y,2)
                        ctrl(k1,:) = bsxfun(@plus, [1 2 3], b(k1).XOffset');
                        ydt(k1,:) = b(k1).YData;
                    end
                    hold on
                    errorbar(ctrl, ydt, errorplus', '.k')

                     % 1/3
                    [p, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[2]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis([ampspost; ampspostP]' ,[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','LFP Amplitude (mV)','')
                    b(1).CData(1,:) = [0 0.4470 0.7410];
                    b(2).CData(1,:) = [0 0.4470 0.7410]/0.7;
                    b(1).CData(2,:) = [0.8500 0.3250 0.0980];
                    b(2).CData(2,:) = [0.8500 0.3250 0.0980]/0.7;
                    b(1).CData(3,:) = [0.4660 0.6740 0.1880];
                    b(2).CData(3,:) = [0.4660 0.6740 0.1880]/0.7;

                    %% Now PE plotting (new way)

                    euclidnorm = sqrt((ampsdev.^2)+(ampsstd.^2)+(ampsctr.^2));
                    euclidnormP = sqrt((ampsdevP.^2)+(ampsstdP.^2)+(ampsctr.^2));

                    ampsdev = ampsdev ./ euclidnorm;
                    ampsctr = ampsctr ./ euclidnorm;
                    ampsstd = ampsstd ./ euclidnorm;
                    ampsdevP = ampsdevP ./ euclidnormP;
                    ampsstdP = ampsstdP ./ euclidnormP;

                    iPE = ampsdev - ampsctr;
                    iRS = ampsctr - ampsstd;
                    iPEp = ampsdevP - ampsctr;
                    iRSp = ampsctr - ampsstdP;

                    iMM = iPE + iRS;
                    iMMp = iPEp + iRSp;

                    nexttile
                    violindata = [iPE', iPEp', iRS', iRSp', iMM', iMMp'];

                    [p, tbl, stats] = kruskalwallis(violindata,[],'off');
                    data = multcompare(stats,'display','off');

                    vp = violinplot(violindata, {'iPE','iPE P','iRS','iRS P','iMM','iMM P'});
                    vpColors = [[255 153 0]*0.85; [255 153 0]; [0 153 204]*0.85; [0 153 204]; [229 0 126]*0.85; [229 0 126]];
                    vpColors = vpColors/255;
                    for bw = 1:length(vp)
                        vp(bw).BoxWidth = 0.05;
                        vp(bw).WhiskerPlot.LineWidth = 2;
                        vp(bw).ViolinColor{1} = vpColors(bw,:);
                    end

                    % 1/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,1:2),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1 2]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,3:4),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[3 4]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,5:6),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','Mismatch index','')

                    %%
                    sgtitle(['Channels = ' num2str(length(iPE))])

                    beautify

                    savepath = 'D:\Jazmin\MultichannelDataTanks\HIP\LFP_Summary';
                    if ~exist(savepath,'dir'), mkdir(savepath); 
                    end

                    print(gcf,'-vector','-depsc',[savepath '\LFP_BP4_0.1-30Hz_Matrix' '.eps'])
                     print(gcf,'-vector','-depsc',[savepath '\LFP_BP4_0.1-30Hz_Matrix' '.pdf'])
                    print(gcf,'-dpng','-r150',[savepath '\LFP_BP4_0.1-30Hz_Matrix' '.png'])

DEV = devM(:,devSamps)*1000;
STD = devM(:,stdSamps)*1000;
STD_Af = devM(:,postSamps)*1000;
CTR = (ctrM)*1000;

SEM_DEV = (std((DEV),"omitnan"))./sqrt(length(DEV));
SEM_STD = (std((STD),"omitnan"))./sqrt(length(STD));
SEM_STD_Af = (std((STD_Af),"omitnan"))./sqrt(length(STD_Af));
SEM_CTR = (std((CTR),"omitnan"))./sqrt(length(CTR));

figure

shadedErrorBar(timevector(devSamps),mean(DEV),SEM_DEV,'lineProps','r')
hold on
shadedErrorBar(timevector(stdSamps),  mean(STD),SEM_STD,'lineProps','b')
shadedErrorBar(timevector(postSamps),  mean(STD_Af),SEM_STD_Af,'lineProps','b')
shadedErrorBar(timevector,mean(CTR),SEM_CTR,'lineProps','g')

figure
SEM_DEV = (std((DEV),"omitnan"))./sqrt(length(DEV));
SEM_STD = (std((STD),"omitnan"))./sqrt(length(STD));

shadedErrorBar(timevector(devSamps),mean(DEV),SEM_DEV,'lineProps','r')
hold on
shadedErrorBar(timevector(devSamps),  mean(STD),SEM_STD,'lineProps','b')

figure
MMN = mean(DEV) - mean(STD);
SEM_MMN = SEM_DEV + SEM_STD;
plot(timevector(devSamps), MMN)
shadedErrorBar(timevector(devSamps),MMN,SEM_MMN,'lineProps','k')

for i = 1:626
    [he(i), pe(i)] = ttest2(DEV(:,i), STD(:,i));
end

figure
plot(timevector(devSamps),he)
[c_pvalues, c_alpha, he, extra] = fdr_BH(pe, 0.01);

figure
plot(timevector(devSamps),c_pvalues)
yline(0.05)
yline(0.01)
yline(0.001)
set(gca, 'YScale', 'log')

DEV_P = devMP(:,devSamps)*1000;
figure
PER = mean(devMP)*1000;
hold on
RAN = mean(devM)*1000;

figure
SEM_DEV_P = (std((devMP)*1000,"omitnan"))./sqrt(length(devMP));
SEM_DEV = (std((devM)*1000,"omitnan"))./sqrt(length(devM));
shadedErrorBar(timevector, PER,SEM_DEV_P)
hold on
shadedErrorBar(timevector, RAN,SEM_DEV)

for i = 1:1879
    [heRP(i), peRP(i)] = ttest2((devMP(:,i))*1000, (devM(:,i))*1000);
end
figure
plot(timevector,heRP)
[c_pvalues, c_alpha, heRP, extra] = fdr_BH(peRP, 0.01);
figure
plot(timevector,c_pvalues)
yline(0.05)
yline(0.01)
yline(0.001)
set(gca, 'YScale', 'log')

figure
shadedErrorBar(timevector(stdSamps), RAN(stdSamps),SEM_DEV(stdSamps),'lineProps','b')
hold on
shadedErrorBar(timevector(devSamps), RAN(devSamps),SEM_DEV(devSamps),'lineProps','r')
shadedErrorBar(timevector(postSamps), RAN(postSamps),SEM_DEV(postSamps),'lineProps','b')
shadedErrorBar(timevector(stdSamps), PER(stdSamps),SEM_DEV_P(stdSamps),'lineProps','--b')
shadedErrorBar(timevector(devSamps), PER(devSamps),SEM_DEV_P(devSamps),'lineProps','--r')
shadedErrorBar(timevector(postSamps), PER(postSamps),SEM_DEV_P(postSamps),'lineProps','--b')


                    close all
                     %% Significance between LFP amplitud of DEV, CTR and STD matrix values
    [h p_iMM] = ttest(iMM) ;
    [h p_iMMP] = ttest(iMMp) ;
    [h p_iRS] = ttest(iRS) ;
    [h p_iRSP] = ttest(iRSp) ;
    [h p_iPE] = ttest(iPE) ;
    [h p_iPEP] = ttest(iPEp) ;

    pValues = [p_iMM; p_iRS; p_iPE];

    m = length(pValues); % Number of tests
    alpha = 0.05; % Standard significance level
    correctedAlpha = alpha / m; % Adjusted significance threshold
    correctedPValues = pValues * m;

    % Display results
    table(pValues, correctedPValues, 'VariableNames', {'pValue', 'Bonfe_Cor_pValue'})

                    [p_s_d, tbl, stats] = kruskalwallis([ampsstd; ampsdev]',[],'off');
                    [p_sP_dP, tbl, stats] = kruskalwallis([ampsstdP; ampsdevP]',[],'off');
                    [p_d_C, tbl, stats] = kruskalwallis([ampsdev; ampsctr]',[],'off');
                    [p_dP_C, tbl, stats] = kruskalwallis([ampsdevP; ampsctr]',[],'off');
                    [p_s_C, tbl, stats] = kruskalwallis([ampsstd; ampsctr]',[],'off');
                    [p_sP_C, tbl, stats] = kruskalwallis([ampsstdP; ampsctr]',[],'off');

 pValues = [p_s_d; p_d_C; p_s_C];
 m = length(pValues); % Number of tests
    alpha = 0.05; % Standard significance level
    correctedAlpha = alpha / m; % Adjusted significance threshold
    correctedPValues = pValues * m;

    % Display results
    table(pValues, correctedPValues, 'VariableNames', {'pValue', 'Bonfe_Cor_pValue'})

iMM_mean = nanmean(iMM)
iMM__std = nanstd(iMM);
err_iMM = iMM__std ./ sqrt(size(iMM,2))

iRS_mean = nanmean(iRS)
iRS__std = nanstd(iRS);
err_iRS = iRS__std ./ sqrt(size(iRS,2))


iPE_mean = nanmean(iPE)
iPE__std = nanstd(iPE);
err_iPE = iPE__std ./ sqrt(size(iPE,2))
                    %%  %% significance

for i = 1:size(mean(devM),2)
    [h(i), p(i)] = ttest2(mean(devM),mean(devMP));
end
figure
plot(timevector,h)
[c_pvalues, c_alpha, h, extra] = fdr_BH(p, 0.05);
figure
plot(timevector,h)
[c_pvalues, c_alpha, h] = fwer_holmbonf(p, 0.05);
figure
plot(timevector,h)
 [c_pvalues, c_alpha, h] = fwer_sidak(p, 0.05);
figure
plot(timevector,h)
[c_pvalues, chi2_2k, h, extra] = mt_fisher(p, 0.05);
figure
plot(timevector,h)


 %% MAtrix analysis of significant reponses
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Matrix created before around line 133
                        % devM_S_DG;
                        % devMP_S_DG;
                        % ctrM_S_DG;
                    %% Added plotting for mean of that layer for each electrode location
                    timevector = (1:length(devM_S_DG))/0.626;
                    timevector = timevector-1000;

                    stdSamps = 1:626;
                    devSamps = 626:1252;
                    postSamps = 1252:length(timevector);

                    figure
                    tiledlayout(1,4)
                    lims = [];
                    nexttile
                    plot(timevector(stdSamps),mean(devM_S_DG(:,stdSamps))*1000,'b')
                    hold on
                    plot(timevector(devSamps),mean(devM_S_DG(:,devSamps))*1000,'r')
                    plot(timevector(postSamps),mean(devM_S_DG(:,postSamps))*1000,'k')

                    plot(timevector(stdSamps),mean(devMP_S_DG(:,stdSamps))*1000,'b--')
                    plot(timevector(devSamps),mean(devMP_S_DG(:,devSamps))*1000,'r--')
                    plot(timevector(postSamps),mean(devMP_S_DG(:,postSamps))*1000,'k--')

                    lims = ylim;
                     ylim([-0.06 0.08])

                    nexttile
                    plot(timevector,mean(ctrM_S_DG)*1000,'g')
                    hold on
                    ylim([-0.06 0.08])
                    %% getting bar data PEAK 2 PEAK
                    ampsstd = peak2peak(devM_S_DG(:,156:406)');
                    ampsstdP = peak2peak(devMP_S_DG(:,156:406)');

                    ampsdev = peak2peak(devM_S_DG(:,783:1033)');
                    ampsdevP = peak2peak(devMP_S_DG(:,783:1033)');

                    ampspost = peak2peak(devM_S_DG(:,1409:1659)');
                    ampspostP = peak2peak(devMP_S_DG(:,1409:1659)');

                    ampsctr = peak2peak(ctrM_S_DG(:,156:406)');

                    nexttile
                    x=categorical({'STD';'DEV';'CTR'});
                    x=reordercats(x,{'STD';'DEV';'CTR'});
                    y=[mean(ampsstd) mean(ampsstdP);  mean(ampsdev) mean(ampsdevP);  mean(ampsctr) 0] *1000;
                    errorplus=([std(ampsstd) std(ampsstdP);  std(ampsdev) std(ampsdevP);  std(ampsctr) 0] *1000 ) / sqrt(length(ampsctr));
                    errorminus=errorplus;
                    b = bar(x, y, 0.8, 'FaceColor' , 'flat');
                    clear ctrl ydt
                    for k1 = 1:size(y,2)
                        ctrl(k1,:) = bsxfun(@plus, [1 2 3], b(k1).XOffset');
                        ydt(k1,:) = b(k1).YData;
                    end
                    hold on
                    errorbar(ctrl, ydt, errorplus', '.k')

                     % 1/3
                    [p, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[2]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis([ampspost; ampspostP]' ,[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','LFP Amplitude (mV)','')
                    b(1).CData(1,:) = [0 0.4470 0.7410];
                    b(2).CData(1,:) = [0 0.4470 0.7410]/0.7;
                    b(1).CData(2,:) = [0.8500 0.3250 0.0980];
                    b(2).CData(2,:) = [0.8500 0.3250 0.0980]/0.7;
                    b(1).CData(3,:) = [0.4660 0.6740 0.1880];
                    b(2).CData(3,:) = [0.4660 0.6740 0.1880]/0.7;

                    %% Now PE plotting (new way)

                    euclidnorm = sqrt((ampsdev.^2)+(ampsstd.^2)+(ampsctr.^2));
                    euclidnormP = sqrt((ampsdevP.^2)+(ampsstdP.^2)+(ampsctr.^2));

                    ampsdev = ampsdev ./ euclidnorm;
                    ampsctr = ampsctr ./ euclidnorm;
                    ampsstd = ampsstd ./ euclidnorm;
                    ampsdevP = ampsdevP ./ euclidnormP;
                    ampsstdP = ampsstdP ./ euclidnormP;

                    iPE = ampsdev - ampsctr;
                    iRS = ampsctr - ampsstd;
                    iPEp = ampsdevP - ampsctr;
                    iRSp = ampsctr - ampsstdP;

                    iMM = iPE + iRS;
                    iMMp = iPEp + iRSp;

                    nexttile
                    violindata = [iPE', iPEp', iRS', iRSp', iMM', iMMp'];

                    [p, tbl, stats] = kruskalwallis(violindata,[],'off');
                    data = multcompare(stats,'display','off');

                    vp = violinplot(violindata, {'iPE','iPE P','iRS','iRS P','iMM','iMM P'});
                    vpColors = [[255 153 0]*0.85; [255 153 0]; [0 153 204]*0.85; [0 153 204]; [229 0 126]*0.85; [229 0 126]];
                    vpColors = vpColors/255;
                    for bw = 1:length(vp)
                        vp(bw).BoxWidth = 0.05;
                        vp(bw).WhiskerPlot.LineWidth = 2;
                        vp(bw).ViolinColor{1} = vpColors(bw,:);
                    end

                    % 1/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,1:2),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1 2]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,3:4),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[3 4]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,5:6),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','Mismatch index DG','')

                    %%
                    sgtitle(['Channels = ' num2str(length(iPE))])

                    beautify

                    savepath = 'D:\Jazmin\MultichannelDataTanks\HIP\LFP_Summary';
                    if ~exist(savepath,'dir'), mkdir(savepath); 
                    end

                    print(gcf,'-vector','-depsc',[savepath '\LFP_BP4_0.1-30Hz_Significant_DG' '.eps'])
                     print(gcf,'-vector','-depsc',[savepath '\LFP_BP4_0.1-30Hz_Significant_DG' '.pdf'])
                    print(gcf,'-dpng','-r150',[savepath '\LFP_BP4_0.1-30Hz_Significant_DG' '.png'])
                    print(gcf,'-dpng','-r150',[savepath '\LFP_BP4_0.1-30Hz_Significant_DG' '.fig'])
                    close all
    %% Significance between LFP amplitud of DEV, CTR and STD matrix values
    [h p_iMM_sig_DG] = ttest(iMM) ;
    [h p_iMMP_sig_DG] = ttest(iMMp) ;
    [h p_iRS_sig_DG] = ttest(iRS);
    [h p_iRSP_sig_DG] = ttest(iRSp) ;
    [h p_iPE_sig_DG] = ttest(iPE); 
    [h p_iPEP_sig_DG] = ttest(iPEp);

    pValues = [p_iMM_sig_DG; p_iRS_sig_DG; p_iPE_sig_DG];
    m = length(pValues); % Number of tests
    alpha = 0.05; % Standard significance level
    correctedAlpha = alpha / m; % Adjusted significance threshold
    correctedPValues = pValues * m;

    % Display results
    table(pValues, correctedPValues, 'VariableNames', {'pValue', 'Bonfe_Cor_pValue'})

                    [p_s_dS_DG, tbl, stats] = kruskalwallis([ampsstd; ampsdev]',[],'off');
                    [p_sP_dPS_DG, tbl, stats] = kruskalwallis([ampsstdP; ampsdevP]',[],'off');
                    [p_d_CS_DG, tbl, stats] = kruskalwallis([ampsdev; ampsctr]',[],'off');
                    [p_dP_CS_DG, tbl, stats] = kruskalwallis([ampsdevP; ampsctr]',[],'off');
                    [p_s_CS_DG, tbl, stats] = kruskalwallis([ampsstd; ampsctr]',[],'off');
                    [p_sP_CS_DG, tbl, stats] = kruskalwallis([ampsstdP; ampsctr]',[],'off');
                    [p_s_sp_DG, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
                    [p_d_dp_DG, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');
        
    pValues = [p_s_dS_DG; p_d_CS_DG; p_s_CS_DG];
    m = length(pValues); % Number of tests
    alpha = 0.05; % Standard significance level
    correctedAlpha = alpha / m; % Adjusted significance threshold
    correctedPValues = pValues * m;

    % Display results
    table(pValues, correctedPValues, 'VariableNames', {'pValue', 'Bonfe_Cor_pValue'})

iMM_mean = nanmean(iMM)
iMM__std = nanstd(iMM);
err_iMM = iMM__std ./ sqrt(size(iMM,2))

iRS_mean = nanmean(iRS)
iRS__std = nanstd(iRS);
err_iRS = iRS__std ./ sqrt(size(iRS,2))


iPE_mean = nanmean(iPE)
iPE__std = nanstd(iPE);
err_iPE = iPE__std ./ sqrt(size(iPE,2))
     %% significance

for i = 1:size(mean(devM_S_DG),2)
    [h(i), p(i)] = ttest2(mean(devM_S_DG),mean(devMP_S_DG));
end
figure
plot(timevector,h)
[c_pvalues, c_alpha, h, extra] = fdr_BH(p, 0.05);
figure
plot(timevector,h)
[c_pvalues, c_alpha, h] = fwer_holmbonf(p, 0.05);
figure
plot(timevector,h)
 [c_pvalues, c_alpha, h] = fwer_sidak(p, 0.05);
figure
plot(timevector,h)
[c_pvalues, chi2_2k, h, extra] = mt_fisher(p, 0.05);
figure
plot(timevector,h)

%%
timevector = (1:length(devM_S_DG))/0.626;
                    timevector = timevector-1000;

                    stdSamps = 1:626;
                    devSamps = 626:1251;
                    postSamps = 1252:length(timevector);

                    figure
                    tiledlayout(1,4)
                    lims = [];
                    nexttile
                    plot(timevector(stdSamps),mean(devM_S_DG(:,stdSamps))*1000,'b')
                    hold on
                    plot(timevector(devSamps),mean(devM_S_DG(:,devSamps))*1000,'r')
                    plot(timevector(postSamps),mean(devM_S_DG(:,postSamps))*1000,'k')

                    plot(timevector(stdSamps),mean(devMP_S_DG(:,stdSamps))*1000,'b--')
                    plot(timevector(devSamps),mean(devMP_S_DG(:,devSamps))*1000,'r--')
                    plot(timevector(postSamps),mean(devMP_S_DG(:,postSamps))*1000,'k--')

                    lims = ylim;
                     ylim([-0.06 0.08])

                    nexttile
                    plot(timevector,mean(ctrM_S_DG)*1000,'g')
                    hold on


DEV_DG = devM_S_DG(:,devSamps)*1000;
STD_DG = devM_S_DG(:,stdSamps)*1000;

sum(sum(isnan(DEV_DG)))
sum(sum(isnan(STD_DG)))

figure
SEM_DEV_DG = (std((DEV_DG),"omitnan"))./sqrt(length(DEV_DG));
SEM_STD_DG = (std((STD_DG),"omitnan"))./sqrt(length(STD_DG));

shadedErrorBar(timevector(devSamps),mean(DEV_DG),SEM_DEV_DG,'lineProps','r')
hold on
shadedErrorBar(timevector(devSamps),  mean(STD_DG),SEM_STD_DG,'lineProps','b')

figure
MMN_DG = mean(DEV_DG) - mean(STD_DG);
SEM_MMN_DG = SEM_DEV_DG + SEM_STD_DG;
plot(timevector(devSamps), MMN_DG)
shadedErrorBar(timevector(devSamps),MMN_DG,SEM_MMN_DG,'lineProps','k')

for i = 1:626
    [he(i), pe(i)] = ttest2(DEV_DG(:,i), STD_DG(:,i));
end

figure
plot(timevector(devSamps),he)
[c_pvalues, c_alpha, he, extra] = fdr_BH(pe, 0.01);

figure
plot(timevector(devSamps),c_pvalues)
yline(0.05)
yline(0.01)
yline(0.001)
set(gca, 'YScale', 'log')
saveas(gcf,'p_R_P.pdf');

DEV_DG_P = devMP_S_DG(:,devSamps)*1000;
figure
PER_DG = mean(devMP_S_DG)*1000;
hold on
RAN_DG = mean(devM_S_DG)*1000;

figure
SEM_DEV_P = (std((devMP_S_DG)*1000,"omitnan"))./sqrt(length(devMP_S_DG));
SEM_DEV = (std((devM_S_DG)*1000,"omitnan"))./sqrt(length(devM_S_DG));
shadedErrorBar(timevector, PER_DG,SEM_DEV_P)
hold on
shadedErrorBar(timevector, RAN_DG,SEM_DEV)

for i = 1:1879
    [heRP(i), peRP(i)] = ttest2((devMP_S_DG(:,i))*1000, (devM_S_DG(:,i))*1000);
end
figure
plot(timevector,heRP)
[c_pvalues, c_alpha, heRP, extra] = fdr_BH(peRP, 0.01);
figure
plot(timevector,c_pvalues)
yline(0.05)
yline(0.01)
yline(0.001)
set(gca, 'YScale', 'log')

figure
shadedErrorBar(timevector(stdSamps), RAN_DG(stdSamps),SEM_DEV(stdSamps),'lineProps','b')
hold on
shadedErrorBar(timevector(devSamps), RAN_DG(devSamps),SEM_DEV(devSamps),'lineProps','r')
shadedErrorBar(timevector(postSamps), RAN_DG(postSamps),SEM_DEV(postSamps),'lineProps','b')
shadedErrorBar(timevector(stdSamps), PER_DG(stdSamps),SEM_DEV_P(stdSamps),'lineProps','--b')
shadedErrorBar(timevector(devSamps), PER_DG(devSamps),SEM_DEV_P(devSamps),'lineProps','--r')
shadedErrorBar(timevector(postSamps), PER_DG(postSamps),SEM_DEV_P(postSamps),'lineProps','--b')

STD_Af_DG = (devM_S_DG(:,postSamps))*1000;
CTR_DG = (ctrM_S_DG)*1000;

SEM_DEV = (std((DEV_DG),"omitnan"))./sqrt(length(DEV_DG));
SEM_STD = (std((STD_DG),"omitnan"))./sqrt(length(STD_DG));
SEM_STD_Af = (std((STD_Af_DG),"omitnan"))./sqrt(length(STD_Af_DG));
SEM_CTR = (std((CTR_DG),"omitnan"))./sqrt(length(CTR_DG));

figure

shadedErrorBar(timevector(devSamps),mean(DEV_DG),SEM_DEV,'lineProps','r')
hold on
shadedErrorBar(timevector(stdSamps),  mean(STD_DG),SEM_STD,'lineProps','b')
shadedErrorBar(timevector(postSamps),  mean(STD_Af_DG),SEM_STD_Af,'lineProps','b')
shadedErrorBar(timevector,mean(CTR_DG),SEM_CTR,'lineProps','g')

%% 
 %% MAtrix analysis of significant reponses
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Matrix created before around line 133
                        % devM_S_CA1;
                        % devMP_S_CA1;
                        % ctrM_S_CA1;
                    %% Added plotting for mean of that layer for each electrode location
                    timevector = (1:length(devM_S_CA1))/0.626;
                    timevector = timevector-1000;

                    stdSamps = 1:626;
                    devSamps = 627:1252;
                    postSamps = 1253:length(timevector);

                    figure
                    tiledlayout(1,4)
                    lims = [];
                    nexttile
                    plot(timevector(stdSamps),mean(devM_S_CA1(:,stdSamps))*1000,'b')
                    hold on
                    plot(timevector(devSamps),mean(devM_S_CA1(:,devSamps))*1000,'r')
                    plot(timevector(postSamps),mean(devM_S_CA1(:,postSamps))*1000,'k')

                    plot(timevector(stdSamps),mean(devMP_S_CA1(:,stdSamps))*1000,'b--')
                    plot(timevector(devSamps),mean(devMP_S_CA1(:,devSamps))*1000,'r--')
                    plot(timevector(postSamps),mean(devMP_S_CA1(:,postSamps))*1000,'k--')

                    lims = ylim;
                     ylim([-0.06 0.08])

                    nexttile
                    plot(timevector,mean(ctrM_S_CA1)*1000,'g')
                    hold on
                    ylim([-0.06 0.08])

                    %% getting bar data PEAK 2 PEAK
                    ampsstd = peak2peak(devM_S_CA1(:,156:406)');
                    ampsstdP = peak2peak(devMP_S_CA1(:,156:406)');

                    ampsdev = peak2peak(devM_S_CA1(:,783:1033)');
                    ampsdevP = peak2peak(devMP_S_CA1(:,783:1033)');

                    ampspost = peak2peak(devM_S_CA1(:,1409:1659)');
                    ampspostP = peak2peak(devMP_S_CA1(:,1409:1659)');

                    ampsctr = peak2peak(ctrM_S_CA1(:,156:406)');

                    nexttile
                    x=categorical({'STD';'DEV';'CTR'});
                    x=reordercats(x,{'STD';'DEV';'CTR'});
                    y=[mean(ampsstd) mean(ampsstdP);  mean(ampsdev) mean(ampsdevP);  mean(ampsctr) 0] *1000;
                    errorplus=([std(ampsstd) std(ampsstdP);  std(ampsdev) std(ampsdevP);  std(ampsctr) 0] *1000 ) / sqrt(length(ampsctr));
                    errorminus=errorplus;
                    b = bar(x, y, 0.8, 'FaceColor' , 'flat');
                    clear ctrl ydt
                    for k1 = 1:size(y,2)
                        ctrl(k1,:) = bsxfun(@plus, [1 2 3], b(k1).XOffset');
                        ydt(k1,:) = b(k1).YData;
                    end
                    hold on
                    errorbar(ctrl, ydt, errorplus', '.k')

                     % 1/3
                    [p, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[2]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis([ampspost; ampspostP]' ,[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','LFP Amplitude (mV)','')
                    b(1).CData(1,:) = [0 0.4470 0.7410];
                    b(2).CData(1,:) = [0 0.4470 0.7410]/0.7;
                    b(1).CData(2,:) = [0.8500 0.3250 0.0980];
                    b(2).CData(2,:) = [0.8500 0.3250 0.0980]/0.7;
                    b(1).CData(3,:) = [0.4660 0.6740 0.1880];
                    b(2).CData(3,:) = [0.4660 0.6740 0.1880]/0.7;

                    %% Now PE plotting (new way)

                    euclidnorm = sqrt((ampsdev.^2)+(ampsstd.^2)+(ampsctr.^2));
                    euclidnormP = sqrt((ampsdevP.^2)+(ampsstdP.^2)+(ampsctr.^2));

                    ampsdev = ampsdev ./ euclidnorm;
                    ampsctr = ampsctr ./ euclidnorm;
                    ampsstd = ampsstd ./ euclidnorm;
                    ampsdevP = ampsdevP ./ euclidnormP;
                    ampsstdP = ampsstdP ./ euclidnormP;

                    iPE = ampsdev - ampsctr;
                    iRS = ampsctr - ampsstd;
                    iPEp = ampsdevP - ampsctr;
                    iRSp = ampsctr - ampsstdP;

                    iMM = iPE + iRS;
                    iMMp = iPEp + iRSp;

                    nexttile
                    violindata = [iPE', iPEp', iRS', iRSp', iMM', iMMp'];

                    [p, tbl, stats] = kruskalwallis(violindata,[],'off');
                    data = multcompare(stats,'display','off');

                    vp = violinplot(violindata, {'iPE','iPE P','iRS','iRS P','iMM','iMM P'});
                    vpColors = [[255 153 0]*0.85; [255 153 0]; [0 153 204]*0.85; [0 153 204]; [229 0 126]*0.85; [229 0 126]];
                    vpColors = vpColors/255;
                    for bw = 1:length(vp)
                        vp(bw).BoxWidth = 0.05;
                        vp(bw).WhiskerPlot.LineWidth = 2;
                        vp(bw).ViolinColor{1} = vpColors(bw,:);
                    end

                    % 1/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,1:2),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[1 2]},[data(1,6)])
                    end
                    % 2/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,3:4),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[3 4]},[data(1,6)])
                    end
                    % 3/3
                    [p, tbl, stats] = kruskalwallis(violindata(:,5:6),[],'off');
                    data = multcompare(stats,'display','off');
                    if data(1,6) < (0.05/3)
                        sigstar({[5 6]},[data(1,6)])
                    end

                    figtexts('','Mismatch index CA1','')

                    %%
                    sgtitle(['Channels = ' num2str(length(iPE))])

                    beautify

                    savepath = 'D:\Jazmin\MultichannelDataTanks\HIP\LFP_Summary';
                    if ~exist(savepath,'dir'), mkdir(savepath); 
                    end

                    print(gcf,'-vector','-depsc',[savepath '\LFP_BP4_0.1-30Hz_Significant_CA1' '.eps'])
                    print(gcf,'-vector','-depsc',[savepath '\LFP_BP4_0.1-30Hz_Significant_CA1' '.pdf'])
                    print(gcf,'-dpng','-r150',[savepath '\LFP_BP4_0.1-30Hz_Significant_CA1' '.png'])
                    print(gcf,'-dpng','-r150',[savepath '\LFP_BP4_0.1-30Hz_Significant_CA1' '.fig'])
                    close all

DEV_CA1 = devM_S_CA1(:,devSamps)*1000;
STD_CA1 = devM_S_CA1(:,stdSamps)*1000;
DEV_CA1_P = devMP_S_CA1(:,devSamps)*1000;

sum(sum(isnan(DEV_CA1)))
sum(sum(isnan(STD_CA1)))
sum(sum(isnan(DEV_CA1_P)))

figure
SEM_DEV_CA1 = (std((DEV_CA1),"omitnan"))./sqrt(length(DEV_CA1));
SEM_STD_CA1 = (std((STD_CA1),"omitnan"))./sqrt(length(STD_CA1));

shadedErrorBar(timevector(devSamps),mean(DEV_CA1),SEM_DEV_CA1,'lineProps','r')
hold on
shadedErrorBar(timevector(devSamps),mean(STD_CA1),SEM_STD_CA1,'lineProps','b')

figure
MMN_CA1 = mean(DEV_CA1) - mean(STD_CA1);
SEM_MMN_CA1 = SEM_DEV_CA1 + SEM_STD_CA1;
plot(timevector(devSamps), MMN_CA1)
shadedErrorBar(timevector(devSamps),MMN_CA1,SEM_MMN_CA1,'lineProps','k')

for i = 1:626
    [he(i), pe(i)] = ttest2(DEV_CA1(:,i), STD_CA1(:,i));
end

% [h, p] = ttest2(DEV_DG, STD_DG);

figure
plot(timevector(devSamps),he)
[c_pvalues, c_alpha, he, extra] = fdr_BH(pe, 0.01);

figure
plot(timevector(devSamps),c_pvalues)
yline(0.05)
yline(0.01)
yline(0.001)
set(gca, 'YScale', 'log')
saveas(gcf,'p_R_P.pdf');

figure
PER_CA1 = mean(devMP_S_CA1)*1000;
hold on
RAN_CA1 = mean(devM_S_CA1)*1000;

figure
SEM_DEV_P = (std((devMP_S_CA1)*1000,"omitnan"))./sqrt(length(devMP_S_CA1));
SEM_DEV = (std((devM_S_CA1)*1000,"omitnan"))./sqrt(length(devM_S_CA1));
shadedErrorBar(timevector, PER_CA1,SEM_DEV_P)
hold on
shadedErrorBar(timevector, RAN_CA1,SEM_DEV)

for i = 1:1879
    [heRP(i), peRP(i)] = ttest2((devMP_S_CA1(:,i))*1000, (devM_S_CA1(:,i))*1000);
end
figure
plot(timevector,heRP)
[c_pvalues, c_alpha, heRP, extra] = fdr_BH(peRP, 0.01);
figure
plot(timevector,c_pvalues)
yline(0.05)
yline(0.01)
yline(0.001)
set(gca, 'YScale', 'log')

figure
shadedErrorBar(timevector(stdSamps), RAN_CA1(stdSamps),SEM_DEV(stdSamps),'lineProps','b')
hold on
shadedErrorBar(timevector(devSamps), RAN_CA1(devSamps),SEM_DEV(devSamps),'lineProps','r')
shadedErrorBar(timevector(postSamps), RAN_CA1(postSamps),SEM_DEV(postSamps),'lineProps','b')
shadedErrorBar(timevector(stdSamps), PER_CA1(stdSamps),SEM_DEV_P(stdSamps),'lineProps','--b')
shadedErrorBar(timevector(devSamps), PER_CA1(devSamps),SEM_DEV_P(devSamps),'lineProps','--r')
shadedErrorBar(timevector(postSamps), PER_CA1(postSamps),SEM_DEV_P(postSamps),'lineProps','--b')

STD_Af_CA1 = (devM_S_CA1(:,postSamps))*1000;
CTR_CA1 = (ctrM_S_CA1)*1000;
                  

SEM_DEV = (std((DEV_CA1),"omitnan"))./sqrt(length(DEV_CA1));
SEM_STD = (std((STD_CA1),"omitnan"))./sqrt(length(STD_CA1));
SEM_STD_Af = (std((STD_Af_CA1),"omitnan"))./sqrt(length(STD_Af_CA1));
SEM_CTR = (std((CTR_CA1),"omitnan"))./sqrt(length(CTR_CA1));

figure

shadedErrorBar(timevector(devSamps),mean(DEV_CA1),SEM_DEV,'lineProps','r')
hold on
shadedErrorBar(timevector(stdSamps),  mean(STD_CA1),SEM_STD,'lineProps','b')
shadedErrorBar(timevector(postSamps),  mean(STD_Af_CA1),SEM_STD_Af,'lineProps','b')
shadedErrorBar(timevector,mean(CTR_CA1),SEM_CTR,'lineProps','g')

    %% Significance between LFP amplitud of DEV, CTR and STD matrix values
    [h p_iMM_sig_CA1] = ttest(iMM) ;
    [h p_iMMP_sig_CA1] = ttest(iMMp) ;
    [h p_iRS_sig_CA1] = ttest(iRS) ;
    [h p_iRSP_sig_CA1] = ttest(iRSp) ;
    [h p_iPE_sig_CA1] = ttest(iPE) ;
    [h p_iPEP_sig_CA1] = ttest(iPEp) ;

    pValues = [p_iMM_sig_CA1; p_iRS_sig_CA1; p_iPE_sig_CA1];
    m = length(pValues); % Number of tests
    alpha = 0.05; % Standard significance level
    correctedAlpha = alpha / m; % Adjusted significance threshold
    correctedPValues = pValues * m;

    % Display results
    table(pValues, correctedPValues, 'VariableNames', {'pValue', 'Bonfe_Cor_pValue'})


                    [p_s_dS_CA1, tbl, stats] = kruskalwallis([ampsstd; ampsdev]',[],'off');
                    [p_sP_dPS_CA1, tbl, stats] = kruskalwallis([ampsstdP; ampsdevP]',[],'off');
                    [p_d_CS_CA1, tbl, stats] = kruskalwallis([ampsdev; ampsctr]',[],'off');
                    [p_dP_CS_CA1, tbl, stats] = kruskalwallis([ampsdevP; ampsctr]',[],'off');
                    [p_s_CS_CA1, tbl, stats] = kruskalwallis([ampsstd; ampsctr]',[],'off');
                    [p_sP_CS_CA1, tbl, stats] = kruskalwallis([ampsstdP; ampsctr]',[],'off');
                    [p_s_sp_CA1, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
                    [p_d_dp_CA1, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');

    pValues = [p_s_dS_CA1; p_d_CS_CA1; p_s_CS_CA1];
    m = length(pValues); % Number of tests
    alpha = 0.05; % Standard significance level
    correctedAlpha = alpha / m; % Adjusted significance threshold
    correctedPValues = pValues * m;

    % Display results
    table(pValues, correctedPValues, 'VariableNames', {'pValue', 'Bonfe_Cor_pValue'})

iMM_mean = nanmean(iMM)
iMM__std = nanstd(iMM);
err_iMM = iMM__std ./ sqrt(size(iMM,2))

iRS_mean = nanmean(iRS)
iRS__std = nanstd(iRS);
err_iRS = iRS__std ./ sqrt(size(iRS,2))


iPE_mean = nanmean(iPE)
iPE__std = nanstd(iPE);
err_iPE = iPE__std ./ sqrt(size(iPE,2))
     %% significance

for i = 1:size(mean(devM_S_CA1),2)
    [h(i), p(i)] = ttest2(mean(devM_S_CA1),mean(devMP_S_CA1));
end
figure
plot(timevector,h)
[c_pvalues, c_alpha, h, extra] = fdr_BH(p, 0.05);
figure
plot(timevector,h)
[c_pvalues, c_alpha, h] = fwer_holmbonf(p, 0.05);
figure
plot(timevector,h)
 [c_pvalues, c_alpha, h] = fwer_sidak(p, 0.05);
figure
plot(timevector,h)
[c_pvalues, chi2_2k, h, extra] = mt_fisher(p, 0.05);
figure
plot(timevector,h)
    %% 
    % Remove electrode locations with very small LFPs
    RemAmps = peak2peak(dev(:,783:1033)');
    RemIdx = RemAmps < 0.00007; % 0.00007 does okay (0.07mv)
    dev(RemIdx,:) = [];
    devP(RemIdx,:) = [];
    ctr(RemIdx,:) = [];

    %% Added plotting for mean of that layer for each electrode location
    timevector = (1:length(dev))/0.626;
    timevector = timevector-1000;

    stdSamps = 1:626;
    devSamps = 626:1252;
    postSamps = 1252:length(timevector);

    figure
    tiledlayout(1,4)
    lims = [];
    nexttile
    plot(timevector(stdSamps),mean(dev(:,stdSamps))*1000,'b')
    hold on
    plot(timevector(devSamps),mean(dev(:,devSamps))*1000,'r')
    plot(timevector(postSamps),mean(dev(:,postSamps))*1000,'k')

    plot(timevector(stdSamps),mean(devP(:,stdSamps))*1000,'b--')
    plot(timevector(devSamps),mean(devP(:,devSamps))*1000,'r--')
    plot(timevector(postSamps),mean(devP(:,postSamps))*1000,'k--')

    lims = ylim;
    ylim([-0.06 0.08])

    plot([0 75],[lims(2)*0.9 lims(2)*0.9],'r')
    plot([0 75]-1000,[lims(2)*0.9 lims(2)*0.9],'b')
    plot([0 75]+1000,[lims(2)*0.9 lims(2)*0.9],'k')
    set(gca,'XLim',[-1000 2000])

    figtexts('ODD','Amplitude (mV)','Time (ms)')

    nexttile
    plot(timevector,mean(ctr)*1000,'g')
    hold on
    ylim([-0.06 0.08])
    figtexts('CTR','','')
    set(gca,'XLim',[-1000 2000])

 %% getting bar data PEAK 2 PEAK
    ampsstd = peak2peak(dev(:,156:406)');
    ampsstdP = peak2peak(devP(:,156:406)');

    ampsdev = peak2peak(dev(:,783:1033)');
    ampsdevP = peak2peak(devP(:,783:1033)');

    ampspost = peak2peak(dev(:,1409:1659)');
    ampspoststdP = peak2peak(devP(:,1409:1659)');

    ampsctr = peak2peak(ctr(:,156:406)');

    nexttile
    x=categorical({'STD';'DEV';'CTR'});
    x=reordercats(x,{'STD';'DEV';'CTR'});
    y=[mean(ampsstd) mean(ampsstdP);  mean(ampsdev) mean(ampsdevP);  mean(ampsctr) 0] *1000;
    errorplus=([std(ampsstd) std(ampsstdP);  std(ampsdev) std(ampsdevP);  std(ampsctr) 0] *1000 ) / sqrt(length(ampsctr));
    errorminus=errorplus;
    b = bar(x, y, 0.8, 'FaceColor' , 'flat');
    clear ctrl ydt
    for k1 = 1:size(y,2)
        ctrl(k1,:) = bsxfun(@plus, [1 2 3], b(k1).XOffset');
        ydt(k1,:) = b(k1).YData;
    end
    hold on
    errorbar(ctrl, ydt, errorplus', '.k')


    % 1/3
    [p, tbl, stats] = kruskalwallis([ampsstd; ampsstdP]',[],'off');
    data = multcompare(stats,'display','off');
    if data(1,6) < (0.05/3)
       sigstar({[1]},[data(1,6)])
    end
    % 2/3
    [p, tbl, stats] = kruskalwallis([ampsdev; ampsdevP]',[],'off');
    data = multcompare(stats,'display','off');

    if data(1,6) < (0.05/3)
       sigstar({[2]},[data(1,6)])
    end

    figtexts('','LFP Amplitude (mV)','')
    b(1).CData(1,:) = [0 0.4470 0.7410];
    b(2).CData(1,:) = [0 0.4470 0.7410];
    b(1).CData(2,:) = [0.8500 0.3250 0.0980];
    b(2).CData(2,:) = [0.8500 0.3250 0.0980];
    b(1).CData(3,:) = [0.4660 0.6740 0.1880];

    %% Now PE plotting (new way)

    euclidnorm = sqrt((ampsdev.^2)+(ampsstd.^2)+(ampsctr.^2));
    euclidnormP= sqrt((ampsdevP.^2)+(ampsstdP.^2)+(ampsctr.^2));

    ampsdev = ampsdev ./ euclidnorm;
    ampsctr = ampsctr ./ euclidnorm;
    ampsstd = ampsstd ./ euclidnorm;
    ampsdevP = ampsdevP ./ euclidnormP;
    ampsctrP = ampsctr./ euclidnormP;
    ampsstdP = ampsstdP ./ euclidnormP;

    iPE = ampsdev - ampsctr;
    iRS = ampsctr - ampsstd;

    iRSp = ampsctrP- ampsstdP;
    iPEp = ampsdevP - ampsctrP;

    iMM = iPE + iRS;
    iMMp = iPEp + iRSp;

    nexttile
    violindata = [iPE', iPEp', iRS', iRSp', iMM', iMMp'];

    vp = violinplot(violindata, {'iPE','iPE P','iRS','iRS P','iMM','iMM P'});

    vpColors = [[255 153 0]*0.85; [255 153 0]; [0 153 204]*0.85; [0 153 204]; [229 0 126]*0.85; [229 0 126]];
    vpColors = vpColors/255;
    for bw = 1:length(vp)
        vp(bw).BoxWidth = 0.05;
        vp(bw).WhiskerPlot.LineWidth = 2;
        vp(bw).ViolinColor{1} = vpColors(bw,:);
    end

    % Stats
    if signrank(violindata(:,1),violindata(:,2)) < (0.05/3)
        sigstar({[1 2]})
    end
    if signrank(violindata(:,3),violindata(:,4)) < (0.05/3)
        sigstar({[3 4]})
    end
    if signrank(violindata(:,5),violindata(:,6)) < (0.05/3)
        sigstar({[5 6]})
    end

    figtexts('','Mismatch index','')

    %%
    sgtitle(['Electrode locations = ' num2str(count)])
    beautify

    savepath = 'D:\Jazmin\MultichannelDataTanks\HIP\LFP_Summary';
    if ~exist(savepath,'dir'), mkdir(savepath); end

    print(gcf,'-vector','-depsc',[savepath '\LFP' '.eps'])
    print(gcf,'-dpng','-r150',[savepath '\LFP' '.png'])
    close all
    disp([layer{la} ' done'])

end


function figtexts(tit,ylab,xlab)
% Just easier to have one line to set title, y and x
% Leave empty string for no text e.g. 
% figtexts('','MYYAXIS','')
title(gca,tit)
ylabel(gca,ylab)
xlabel(gca,xlab)
end