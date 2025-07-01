% Analyse Single unit oddball data from sorted bst data
% Adam Hockley 22/11/2023
% Modified by Jazmin Sanchez 2024


clear all
clc

Main_path = 'D:\Jazmin\MultichannelDataTanks\HIP\';
tanks = {'24_023','24_024','24_030'};
Sorter = 'KS4';
warning('off','MATLAB:table:RowsAddedExistingVars')
SR = 1000; % SOA
SpikeWindow = []; % Only take spike from a set time window (ms)

spiketable = table({},[],{},{},{},{},{},{},{},[],'VariableNames',["Block","Unit","Spikes","SpikesPdev","SpikesMS","SpikesCASCASC","SpikesCASCDES","oddData","Shapes","SFR"]);

for ta = 1:length(tanks)

    % Get electrode positions
    Positions = allfolders([Main_path tanks{ta}]);

    for pos = 1:length(Positions)

        tank_path = [Main_path tanks{ta} '\' Positions{pos}];

        % Use this bit of code to get all Stim blocks in the tank
        Blocks = allfolders(tank_path);
        Blocks = Blocks(contains(Blocks,'STIM'));
        Blocks = Blocks(~contains(Blocks,'1000'));

        for bl = 1:length(Blocks)

            block_path = [tank_path '\' Blocks{bl}];

            if exist([block_path '\bst_' Sorter '.mat'])
                load([block_path '\bst_' Sorter '.mat']) % load the bst

                Units = unique(bst.Spikes.unit);

                for un = 1:length(Units)

                    %% Oddball plotting
                    oddData = BST_ODD3_Jazmin(bst,Units(un));

                    %% Get spont rate
                    SelectedTrials = 1:20; % first 20 trials are silent
                    SpikeTimes = SR*BST_GS3(bst,SelectedTrials,Units(un));
                    spontratehz = length(SpikeTimes) / bst.Epocs.TSOff.wfrq(20);

                    %% Get histogram spike windows etc
                    if isempty(SpikeWindow)
                        SpikeWindow = [0 round((bst.Epocs.TSOn.wfrq(20) - bst.Epocs.TSOn.wfrq(19))*SR)*3]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hard coded sampling rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                    histwidth = 25; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histogram width %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    histedges = -SR:histwidth:SR;
                    histcenters = histedges(1:end-1)+histwidth/2;
                    histcentersinterp = histcenters; % get new xaxis for interpolated data
                    %                 histcentersinterp = linspace(histcenters(1),histcenters(end),1000); % get new xaxis for interpolated data

                    %%
                    DEVa = vertcat(oddData.('ODD-ASC-DEV'){:});
                    DEVb = vertcat(oddData.('ODD-DES-DEV'){:});
                    spikes = [DEVa; DEVb];
                    spikes = histcounts(spikes,histedges);
                    spikes = spikes/((histwidth/SR)*80); % convert to Hz

                    DEVa = vertcat(oddData.('ODD-ASC-DEV-P'){:});
                    DEVb = vertcat(oddData.('ODD-DES-DEV-P'){:});
                    spikesPdev = [DEVa;DEVb];
                    spikesPdev = histcounts(spikesPdev,histedges);
                    spikesPdev = spikesPdev/((histwidth/SR)*80); % convert to Hz

                    MSa = vertcat(oddData.('MS-F1'){:});
                    MSb = vertcat(oddData.('MS-F2'){:});
                    spikesMS = [MSa;MSb];
                    spikesMS = histcounts(spikesMS,histedges);
                    spikesMS = spikesMS/((histwidth/SR)*80); % convert to Hz

                    CASCASCa = vertcat(oddData.('CASC-ASC-F1'){:});
                    CASCASCb = vertcat(oddData.('CASC-ASC-F2'){:});
                    spikesCASCASC = [ CASCASCa;CASCASCb];
                    spikesCASCASC = histcounts(spikesCASCASC,histedges);
                    spikesCASCASC = spikesCASCASC/((histwidth/SR)*80); % convert to Hz

                    CASCDESa = vertcat(oddData.('CASC-DES-F1'){:});
                    CASCDESb = vertcat(oddData.('CASC-DES-F2'){:});
                    spikesCASCDES = [ CASCDESa; CASCDESb ];
                    spikesCASCDES = histcounts(spikesCASCDES,histedges);
                    spikesCASCDES = spikesCASCDES/((histwidth/SR)*80); % convert to Hz

                    %% Add to data table
                    textloc = strfind(block_path,'HIP');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    spiketable.Block{end+1} = block_path(textloc+4:end);
                    spiketable.Unit(end) = Units(un);
                    spiketable.Spikes{end} = spikes;
                    spiketable.SpikesPdev{end} = spikesPdev;
                    spiketable.SpikesMS{end} = spikesMS;
                    spiketable.SpikesCASCASC{end} = spikesCASCASC;
                    spiketable.SpikesCASCDES{end} = spikesCASCDES;
                    spiketable.oddData{end} = oddData;
                    spiketable.SFR(end) = spontratehz;

                    %% Plot

                    tiledlayout(1,5);
                    set(gcf,'Visible','off')
                    nexttile([1 2])
                    plot(histcentersinterp,spikes,'k')
                    hold on
                    plot(histcentersinterp,spikesPdev,'k--')
                    xticks(histedges(1:20:end))
                    temp = get(gca,'ylim');
                    ydat(1) = temp(2);
                    figtexts('ODD','Spike rate (Hz)','Time (ms)')

                    nexttile
                    plot(histcentersinterp,spikesMS,'k')
                    hold on
                    xticks(histedges(1:10:end))
                    figtexts('MS','Spike rate (Hz)','Time (ms)')
                    temp = get(gca,'ylim');
                    ydat(2) = temp(2);
                    set(gca,'XLim',[-5000 5000])

                    nexttile
                    plot(histcentersinterp,spikesCASCASC,'k')
                    hold on
                    xticks(histedges(1:10:end))
                    figtexts('CASCASC','Spike rate (Hz)','Time (ms)')
                    temp = get(gca,'ylim');
                    ydat(2) = temp(2);
                    set(gca,'XLim',[-5000 5000])

                    nexttile
                    plot(histcentersinterp,spikesCASCDES,'k')
                    hold on
                    xticks(histedges(1:10:end))
                    figtexts('CASCDESC','Spike rate (Hz)','Time (ms)')
                    temp = get(gca,'ylim');
                    ydat(2) = temp(2);
                    set(gca,'XLim',[-5000 5000])

                    nexttile(1)
                    ylim([0 max(ydat)])
                    plot([-2000 -1925],[max(ydat) max(ydat)]*0.8,'b','linewidth',2)
                    plot([-1000 -925],[max(ydat) max(ydat)]*0.8,'b','linewidth',2)
                    plot([0 75],[max(ydat) max(ydat)]*0.8,'r','linewidth',2)
                    plot([1000 1075],[max(ydat) max(ydat)]*0.8,'k','linewidth',2)
                    nexttile(3)
                    ylim([0 max(ydat)])
                    plot([0 75],[max(ydat) max(ydat)]*0.8,'g','linewidth',2)
                    nexttile(4)
                    ylim([0 max(ydat)])
                    plot([0 75],[max(ydat) max(ydat)]*0.8,'g','linewidth',2)
                    nexttile(5)
                    ylim([0 max(ydat)])
                    plot([0 75],[max(ydat) max(ydat)]*0.8,'g','linewidth',2)

                    beautify

                    %%
                    savpath = [tank_path '\' Blocks{bl} '\PSTH_plots'];
                    if ~exist(savpath,'dir'), mkdir(savpath); end

                    print(gcf,'-vector','-depsc',[savpath '\' num2str(Units(un)) '_PSTH' Sorter '.eps'])
                    print(gcf,'-dpng','-r150',[savpath '\' num2str(Units(un)) '_PSTH' Sorter '.png'])
                    close all

                end
            end
        end
        if exist('textloc')
            disp([block_path(textloc:end) ' done'])
        end
    end
end

% Add unit depth to spike table

for i = 1:height(spiketable)
    Block = spiketable.Block{i};
    slashlocs = strfind(Block,'\');
    load([Main_path Block(1:slashlocs(2)) '\ClusterGoodLocs.mat'])
    unitloc = find(ClusterGood(1,:) == spiketable.Unit(i));
    spiketable.unitloc(i) = double(ClusterGood(2,unitloc));
end

Main_path = 'D:\Jazmin\MultichannelDataTanks\HIP\Marzo2025\';
save([Main_path 'spiketable' Sorter '.mat'],'spiketable') % save the spiketable
save([Main_path 'spiketable_' Sorter '.mat'],'spiketable') % save the spiketable


%%
function data = BST_ODD3(bst,un)
% Creates the oddball table from a bst, which spikes were in which
% conditions and trials etc.

Epoks = bst.Epocs.Values;

% Get the list of frequencies and levels used
SelectedTrials = BST_TS3(bst,'type','ODD-DES-DEV');
F1 = unique(Epoks.wfrq(SelectedTrials(1)));
SelectedTrials = BST_TS3(bst,'type','ODD-ASC-DEV');
F2 = unique(Epoks.wfrq(SelectedTrials(1)));

data = table;
conditions = uniqueCellVector(Epoks.type);
for co = 2:length(conditions)

    SelectedTrials = BST_TS3(bst,'type',conditions{co});
    for i = 1:length(SelectedTrials)

        
        SpikeTimes = 1000 * BST_GS3(bst,SelectedTrials(i),un);
        data.(conditions{co})(i) = {SpikeTimes};
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

%%
function figtexts(tit,ylab,xlab)
% Just easier to have one line to set title, y and x
% Leave empty string for no text e.g. 
% figtexts('','MYYAXIS','')

title(gca,tit)
ylabel(gca,ylab)
xlabel(gca,xlab)
% box off

end