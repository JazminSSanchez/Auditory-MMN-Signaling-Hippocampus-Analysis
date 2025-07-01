function [HistogramsSpikes] = SignificanceWindow_Spikes(spiketable,x1,x2)
%   Detailed explanation goes here

    histwidth = 30; %Width of the histogram bins

    %Bin corresponding to the start of the window of analysis
    x1 = x1./histwidth + 80;
    x1 = round(x1);
    %Bin corresponding to the end of the window of analysis
    x2 = x2./histwidth + 81;
    x2 = round(x2);
    %Bin corresponding to the end of the basal activity window    
    x1b = 100/histwidth + 81;
    x1b = round(x1b);
    %It stores the spikes corresponding to each of the specified labels in
    %the corresponding variable 

    dev2spikes = vertcat(spiketable.('ODD-ASC-DEV'){:});
    dev1spikes = vertcat(spiketable.('ODD-DES-DEV'){:});
    std1spikes = vertcat(spiketable.('ODD-ASC-STD'){:});
    std2spikes = vertcat(spiketable.('ODD-DES-STD'){:});
    dev2spikesP = vertcat(spiketable.('ODD-ASC-DEV-P'){:});
    dev1spikesP = vertcat(spiketable.('ODD-DES-DEV-P'){:});
    std1spikesP = vertcat(spiketable.('ODD-ASC-STD-P'){:});
    std2spikesP = vertcat(spiketable.('ODD-DES-STD-P'){:});
    ctr1caspikes = vertcat(spiketable.('CASC-ASC-F1'){:});
    ctr1cdspikes = vertcat(spiketable.('CASC-DES-F1'){:});
    ctr2caspikes = vertcat(spiketable.('CASC-ASC-F2'){:});
    ctr2cdspikes = vertcat(spiketable.('CASC-DES-F2'){:});
    ctr1MSspikes = vertcat(spiketable.('MS-F1'){:});
    ctr2MSspikes = vertcat(spiketable.('MS-F2'){:});
%% 
    %It generates the bin center values to use the fuction hist
    % histedges = 0:histwidth:1000;
    % histcenters = histedges(1:end-1)+histwidth/2;

    % DEVa = vertcat(oddData.('ODD-ASC-DEV'){:});
    % DEVb = vertcat(oddData.('ODD-DES-DEV'){:});
    % spikes = [DEVa; DEVb];

    histwidth = 25;
    histedges = -2000:histwidth:3000;
    histcenters = histedges(1:end-1)+histwidth/2;
    histcentersinterp = histcenters; % get new xaxis for interpolated data
    
    % dev2spikes = histcounts(dev2spikes,histedges);
    % dev2spikes = dev2spikes/((histwidth/1000)*80);

    tiledlayout(1,14); % Creates a tiled chart layout for displaying multiple plots
    nexttile

    spikeCountstd1 = hist(std1spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatstd1 = 1000 * spikeCountstd1 / (histwidth*(size(std1spikes,1)));
    %Graph
    bar(histcenters,histdatstd1,'k','BarWidth',1);

%     histdatstd1 = histcounts(std1spikes,histedges);
%     histdatstd1 = (histdatstd1/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd1,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD1');
    temp = get(gca,'ylim');
    ydat(1) = temp(2);

    nexttile
    
    spikeCountstd2 = hist(std2spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatstd2 = 1000 * spikeCountstd2 / (histwidth*(size(std2spikes,1)));

    bar(histcenters,histdatstd2,'k','BarWidth',1);

%     histdatstd2 = histcounts(std2spikes,histedges);
%     histdatstd2 = (histdatstd2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD2');
    temp = get(gca,'ylim');
    ydat(2) = temp(2);

    nexttile

    spikeCountctr1 = hist(ctr1caspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr1ca = 1000 * spikeCountctr1 / (histwidth*(size(ctr1caspikes,1)));

    bar(histcenters,histdatctr1ca,'k','BarWidth',1);
    
%     histdatctr1ca = histcounts(ctr1caspikes,histedges);
%     histdatctr1ca = (histdatctr1ca/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr1ca,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR1');
    temp = get(gca,'ylim');
    ydat(3) = temp(2);

    nexttile

    spikeCountctr2 = hist(ctr2caspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr2ca = 1000 * spikeCountctr2 / (histwidth*(size(ctr2caspikes,1)));

    bar(histcenters,histdatctr2ca,'k','BarWidth',1);    

%     histdatctr2ca = histcounts(ctr2caspikes,histedges);
%     histdatctr2ca = (histdatctr2ca/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2ca,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR2');
    temp = get(gca,'ylim');
    ydat(4) = temp(2);

     nexttile

    spikeCountctr1cd = hist(ctr1cdspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr1cd = 1000 * spikeCountctr1cd / (histwidth*(size(ctr1cdspikes,1)));

    bar(histcenters,histdatctr1cd,'k','BarWidth',1);

%     histdatctr1cd = histcounts(ctr1cdspikes,histedges);
%     histdatctr1cd = (histdatctr1cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr1cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR1');
    temp = get(gca,'ylim');
    ydat(5) = temp(2);

    nexttile

    spikeCountctr2cd = hist(ctr2cdspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr2cd = 1000 * spikeCountctr2cd / (histwidth*(size(ctr2cdspikes,1)));

    bar(histcenters,histdatctr2cd,'k','BarWidth',1);    

%     histdatctr2cd = histcounts(ctr2cdspikes,histedges);
%     histdatctr2cd = (histdatctr2cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR2');
    temp = get(gca,'ylim');
    ydat(6) = temp(2);

    nexttile

    spikeCountctr1MS = hist(ctr1MSspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr1MS = 1000 * spikeCountctr1MS / (histwidth*(size(ctr1MSspikes,1)));

    bar(histcenters,histdatctr1MS,'k','BarWidth',1);    

%     histdatctr2cd = histcounts(ctr2cdspikes,histedges);
%     histdatctr2cd = (histdatctr2cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR1 MS');
    temp = get(gca,'ylim');
    ydat(7) = temp(2);
        
    nexttile

    spikeCountctr2MS = hist(ctr2MSspikes,histcenters);
    % normalize and convert to density in spk/s
    histdatctr2MS = 1000 * spikeCountctr2MS / (histwidth*(size(ctr2MSspikes,1)));

    bar(histcenters,histdatctr2MS,'k','BarWidth',1);    

%     histdatctr2cd = histcounts(ctr2cdspikes,histedges);
%     histdatctr2cd = (histdatctr2cd/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatctr2cd,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('CTR2 MS');
    temp = get(gca,'ylim');
    ydat(8) = temp(2);


    nexttile

    spikeCountdev1 = hist(dev1spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatdev1 = 1000 * spikeCountdev1 / (histwidth*(size(dev1spikes,1)));

    bar(histcenters,histdatdev1,'k','BarWidth',1);

%     histdatdev1 = histcounts(dev1spikes,histedges);
%     histdatdev1 = (histdatdev1/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatdev1,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV1');
    temp = get(gca,'ylim');
    ydat(9) = temp(2);

    nexttile

    spikeCountdev2 = hist(dev2spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatdev2 = 1000 * spikeCountdev2 / (histwidth*(size(dev2spikes,1)));

    bar(histcenters,histdatdev2,'k','BarWidth',1);
    
    hold on
    xticks(histedges(1:5:end));
    title('DEV2');
    temp = get(gca,'ylim');
    ydat(10) = temp(2);

    nexttile

    spikeCountstd1P = hist(std1spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatstd1P = 1000 * spikeCountstd1P / (histwidth*(size(std1spikesP,1)));
    %Graph
    bar(histcenters,histdatstd1P,'k','BarWidth',1);

%     histdatstd1 = histcounts(std1spikes,histedges);
%     histdatstd1 = (histdatstd1/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd1,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD1 P');
    temp = get(gca,'ylim');
    ydat(11) = temp(2);

    nexttile
    
    spikeCountstd2P = hist(std2spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatstd2P = 1000 * spikeCountstd2P / (histwidth*(size(std2spikesP,1)));

    bar(histcenters,histdatstd2P,'k','BarWidth',1);

%     histdatstd2 = histcounts(std2spikes,histedges);
%     histdatstd2 = (histdatstd2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatstd2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('STD2 P');
    temp = get(gca,'ylim');
    ydat(12) = temp(2);

    nexttile
    
    spikeCountdev1P = hist(dev1spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatdev1P = 1000 * spikeCountdev1P / (histwidth*(size(dev1spikesP,1)));

    bar(histcenters,histdatdev1P,'k','BarWidth',1);

%     histdatdev2 = histcounts(dev2spikes,histedges);
%     histdatdev2 = (histdatdev2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatdev2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV1 P');
    temp = get(gca,'ylim');
    ydat(13) = temp(2);

    nexttile

    spikeCountdev2P = hist(dev2spikesP,histcenters);
    % normalize and convert to density in spk/s
    histdatdev2P = 1000 * spikeCountdev2P / (histwidth*(size(dev2spikesP,1)));

    bar(histcenters,histdatdev2P,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV2P');
    temp = get(gca,'ylim');
    ydat(14) = temp(2);

    nexttile(1)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(2)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(3)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(4)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(5)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(6)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(7)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(8)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(9)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(10)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(11)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(12)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(13)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    nexttile(14)
    ylim([0 max(ydat)])
    plot([0 75],[max(ydat) max(ydat)]*0.8,'r')
    beautify
% histcenters = table2array(histcenters);
% histdatstd1 = table2array(histdatstd1);
% histdatstd2 = table2array(histdatstd2);
% histdatdev1 = table2array(histdatdev1);
% histdatdev2 = table2array(histdatdev2);
% histdatctr1ca = table2array( histdatctr1ca);
% histdatctr2ca = table2array(histdatctr2ca);
% histdatctr1cd = table2array(histdatctr1cd);
% histdatctr2cd = table2array(histdatctr2cd);
% histdatctr1MS = table2array(histdatctr1MS);
% spikeCountctr2MS = table2array(spikeCountctr2MS);
% histdatstd1P = table2array(histdatstd1P);
% histdatstd2P = table2array(histdatstd2P);
% histdatdev1P = table2array(histdatdev1P);
% histdatdev2P = table2array(histdatdev2P);

HistogramsSpikes = table(histcenters, histdatstd1, histdatstd2, histdatdev1, histdatdev2, histdatctr1ca, histdatctr2ca, histdatctr1cd, histdatctr2cd,histdatctr1MS,spikeCountctr2MS,histdatstd1P,histdatstd2P,histdatdev1P,histdatdev2P);
