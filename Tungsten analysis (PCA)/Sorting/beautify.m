function beautify
% Defaults
set(0,'defaultfigurecolor',[1 1 1]) % white background
set(0,'defaultaxesfontname','arial') % beautify the axes a bit

% General
set(gcf, 'Color', [1,1,1], 'Position', [150 150 600 500]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)

%% Axis specific things, loop through tiles if it's tiled layout
if strcmp(get(gcf().Children, 'Type'), 'tiledlayout')
    grid = gcf().Children.GridSize;
    for i = 1:(grid(1) * grid(2))
        nexttile(i)
        set(gca,'gridcolor',[1 1 1],'linewidth',1) % set the axis color
        box off
        if isnumeric(xticklabels) && length(xticks)>5
            ticks = xticks;
            set(gca, 'XTick', ticks(1:2:end))
        end
        if isnumeric(yticklabels) && length(yticks)>5
            ticks = yticks;
            set(gca, 'YTick', ticks(1:2:end))
        end
        set(gca,'TickDir','out');
        temp = gca;
        temp.Title.FontSize = 14;
        temp.YLabel.FontSize = 14;
        temp.XLabel.FontSize = 14;

    end

    if (grid(2)/4) > grid(1) 
        set(gcf, 'Color', [1,1,1], 'Position', [100 100 1400 500]);
    elseif (grid(2)/2) > grid(1) 
        set(gcf, 'Color', [1,1,1], 'Position', [100 100 1000 500]);
    elseif (grid(1)/4) > grid(2) 
        set(gcf, 'Color', [1,1,1], 'Position', [100 100 500 1400]);
    elseif (grid(1)/2) > grid(2) 
        set(gcf, 'Color', [1,1,1], 'Position', [100 100 500 1000]);
    elseif grid(1)+grid(2) > 4
        set(gcf, 'Color', [1,1,1], 'Position', [100 100 1000 600]);
    end

else
    set(gca,'gridcolor',[1 1 1],'linewidth',1) % set the axis color
    box off
    if isnumeric(xticks)
        ticks = xticks;
        set(gca, 'XTick', ticks(1:2:end))
    end
    if isnumeric(yticks)
        ticks = yticks;
        set(gca, 'YTick', ticks(1:2:end))
    end
    set(gca,'TickDir','out');
    temp = gca;
    temp.Title.FontSize = 14;
    temp.YLabel.FontSize = 14;
    temp.XLabel.FontSize = 14;

end


end

