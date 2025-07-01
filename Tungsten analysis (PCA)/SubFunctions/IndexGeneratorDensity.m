function [NormalizedIndexes histcenters histdatstd1 histdatstd2 histdatdev1 histdatdev2 histdatctr1ca histdatctr2ca histdatctr1cd histdatctr2cd] = IndexGeneratorDensity(spikedata,y1,y2)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here


    y1 = y1./25 + 1;
    y1 = round(y1);
    y2 = y2./25 + 1;
    y2 = round(y2);
    y1b = 100/25 + 1;
    y1b = round(y1b);

    dev2spikes = [spikedata.('ODD-ASC-DEV'){1}];
    dev1spikes = [spikedata.('ODD-DES-DEV'){1}];
    std1spikes = [spikedata.('ODD-ASC-STD'){1}];
    std2spikes = [spikedata.('ODD-DES-STD'){1}];
    ctr1caspikes = [spikedata.('CASC-ASC-F1'){1}];
    ctr1cdspikes = [spikedata.('CASC-DES-F1'){1}];
    ctr2caspikes = [spikedata.('CASC-ASC-F2'){1}];
    ctr2cdspikes = [spikedata.('CASC-DES-F2'){1}];

%% 
    histwidth = 25;
    histedges = 0:histwidth:1000;
    histcenters = histedges(1:end-1)+histwidth/2;

    tiledlayout(1,8);
    nexttile

    spikeCountstd1 = hist(std1spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatstd1 = 1000 * spikeCountstd1 / (histwidth*(size(std1spikes,1)));

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
    ydat(7) = temp(2);

    nexttile

    spikeCountdev2 = hist(dev2spikes,histcenters);
    % normalize and convert to density in spk/s
    histdatdev2 = 1000 * spikeCountdev2 / (histwidth*(size(dev2spikes,1)));

    bar(histcenters,histdatdev2,'k','BarWidth',1);

%     histdatdev2 = histcounts(dev2spikes,histedges);
%     histdatdev2 = (histdatdev2/40)/(histwidth/40); % normalise to Hz
%     bar(histcenters,histdatdev2,'k','BarWidth',1);

    hold on
    xticks(histedges(1:5:end));
    title('DEV2');
    temp = get(gca,'ylim');
    ydat(8) = temp(2);

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

    beautify
%% 
%BasalActivity reduction
    nd1 = size(histdatdev1,2);
    BasalActivityMeanyDev1 = mean(histdatdev1(1:y1b));
    
    ns1 = size(histdatstd1,2);
    BasalActivityMeanySTD1 = mean(histdatstd1(1:y1b));
    
    nd2 = size(histdatdev2,2);
    BasalActivityMeanyDev2 = mean(histdatdev2(1:y1b));
    
    ns2 = size(histdatstd2,2);
    BasalActivityMeanySTD2 = mean(histdatstd2(1:y1b));
    
    nca1 = size(histdatctr1ca,2);
    BasalActivityctr1cas = mean(histdatctr1ca(1:y1b));
    
    ncd1 = size(histdatctr1cd,2);
    BasalActivityMeanyctr1cd = mean(histdatctr1cd(1:y1b));
    
    nca2 = size(histdatctr2ca,2);
    BasalActivityMeanyctr2cas = mean(histdatctr2ca(1:y1b));
    
    ncd2 = size(histdatctr2cd,2);
    BasalActivityMeanyctr2cd = mean(histdatctr2cd(1:y1b));
    
    for i = 1:nd1
        if histdatdev1(i) < BasalActivityMeanyDev1
           histdatdev1(i) = 0;
        else
           histdatdev1(i) =  histdatdev1(i);
        end
    end
    
    for i = 1:nd2
        if histdatdev2(i) < BasalActivityMeanyDev2
           histdatdev2(i) = 0;
        else
           histdatdev2(i) = histdatdev2(i);
        end
    end
    
    for i = 1:ns1
        if histdatstd1(i) < BasalActivityMeanySTD1
           histdatstd1(i) = 0;
        else
           histdatstd1(i) = histdatstd1(i);
        end
    end
    
    for i = 1:ns2
        if histdatstd2(i) < BasalActivityMeanySTD2
           histdatstd2(i) = 0;
        else
            histdatstd2(i) = histdatstd2(i);
        end
    end
    
    for i = 1:nca1
        if histdatctr1ca(i) < BasalActivityctr1cas
           histdatctr1ca(i) = 0;
        else
           histdatctr1ca(i) =  histdatctr1ca(i);
        end
    end
    
    for i = 1:ncd1
        if histdatctr1cd(i) < BasalActivityMeanyctr1cd
           histdatctr1cd(i) = 0;
        else
           histdatctr1cd(i) = histdatctr1cd(i);
        end
    end
    
    for i = 1:nca2
        if histdatctr2ca(i) < BasalActivityMeanyctr2cas
           histdatctr2ca(i) = 0;
        else
           histdatctr2ca(i) = histdatctr2ca(i);
        end
    end
    
    for i = 1:ncd2
        if histdatctr2cd(i) < BasalActivityMeanyctr2cd
           histdatctr2cd(i) = 0;
        else
            histdatctr2cd(i) = histdatctr2cd(i);
        end
    end
    
    % if max(dev1spikes) > max(dev2spikes)
    %     ymax = max(dev1spikes) + 0.3;
    % else
    %     ymax = max(dev2spikes) + 0.3;
    % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

DEV1 = mean(histdatdev1(y1:y2));
STD1 = mean(histdatstd1(y1:y2));

DEV2 = mean(histdatdev2(y1:y2));
STD2 = mean(histdatstd2(y1:y2));

CTRCA1 = mean(histdatctr1ca(y1:y2));
CTRCA2 = mean(histdatctr2ca(y1:y2));

CTRCD1 = mean(histdatctr1cd(y1:y2));
CTRCD2 = mean(histdatctr2cd(y1:y2));

%%%%%%%%%%%%%%%%%%%%%% SI & CSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SI1 = (DEV1-STD1)/(DEV1+STD1);
SI2 = (DEV2-STD2)/(DEV2+STD2);

CSI = ((DEV1+DEV2)-(STD1+STD2))/(DEV1+DEV2+STD1+STD2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nca1 = sqrt( DEV1^2 + STD1^2 + CTRCA1^2);
Nca2 = sqrt( DEV2^2 + STD2^2 + CTRCA2^2);

Ncd1 = sqrt( DEV1^2 + STD1^2 + CTRCD1^2);
Ncd2 = sqrt( DEV2^2 + STD2^2 + CTRCD2^2);

NormDEVca1 = DEV1/Nca1;
NormSTDca1 = STD1/Nca1;

NormDEVcd1 = DEV1/Ncd1;
NormSTDcd1 = STD1/Ncd1;

NormCTRca1 = CTRCA1/Nca1;
NormCTRcd1 = CTRCD1/Ncd1;

NormDEVca2 = DEV2/Nca2;
NormSTDca2 = STD2/Nca2;

NormDEVcd2 = DEV2/Ncd2;
NormSTDcd2 = STD2/Ncd2;

NormCTRca2 = CTRCA2/Nca2;
NormCTRcd2 = CTRCD2/Ncd2;


%%%%%%%%%%%%%%% CASCASC %%%%%%%%%%%%%%%%%

    iMM1ca = NormDEVca1 - NormSTDca1;
    iRS1ca = NormCTRca1 - NormSTDca1;
    iPE1ca = NormDEVca1 - NormCTRca1;


%%%%%%%%%%%%%%% CASCDESC %%%%%%%%%%%%%%%%%

    iMM2cd = NormDEVcd2 - NormSTDcd2;
    iRS2cd = NormCTRcd2 - NormSTDcd2;
    iPE2cd = NormDEVcd2 - NormCTRcd2;


NormalizedIndexes = table(CSI,SI1, iMM1ca, iRS1ca, iPE1ca,SI2, iMM2cd, iRS2cd, iPE2cd, NormDEVca1, NormSTDca1, NormCTRca1, NormDEVcd2, NormSTDcd2, NormCTRcd2, DEV1, STD1, CTRCA1, DEV2, STD2, CTRCD2);

end