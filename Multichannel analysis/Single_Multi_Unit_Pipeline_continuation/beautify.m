function beautify

% General
set(gcf, 'Position', [150 150 600 500]);
set(findall(gcf,'-property','FontSize'),'FontSize',8)
% set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
fontname(gcf,'arial')
set(0,'defaultaxesfontname','arial') % beautify the axes a bit
% set(0,'defaultfigurecolor',[1 1 1]) % white background

%% Axis specific things, loop through tiles if it's tiled layout
if strcmp(get(gcf().Children, 'Type'), 'tiledlayout')

    tlo = get(gcf,'Children');
    tiles = tilenum(tlo.Children);
    %     [row,col] = tilerowcol((tlo.Children));
    grid = gcf().Children.GridSize;

    for i = 1:grid(1)
        for ii = 1:grid(2)
            temploc = tilenum(tlo,i,ii);
            while ~ismember(temploc,tiles)
                temploc = temploc-1;
            end
            tilemat(i,ii) = temploc;
        end
    end

    tlo.TileSpacing = 'compact';
    tlo.Padding = 'compact';

    if max(size(tilemat)) >= 4
        set(gcf, 'Position', [100 100 250*width(tilemat) 200*height(tilemat)]);
    else
        set(gcf, 'Position', [100 100 400*width(tilemat) 350*height(tilemat)]);
    end

    %     if (width(tilemat)/4) > height(tilemat)
    %         set(gcf, 'Color', [1,1,1], 'Position', [100 100 1400 500]);
    %     elseif (width(tilemat)/2) > height(tilemat)
    %         set(gcf, 'Color', [1,1,1], 'Position', [100 100 1000 500]);
    %     elseif (height(tilemat)/4) > width(tilemat)
    %         set(gcf, 'Color', [1,1,1], 'Position', [100 100 500 1400]);
    %     elseif (height(tilemat)/2) > width(tilemat)
    %         set(gcf, 'Color', [1,1,1], 'Position', [100 100 500 1000]);
    %     elseif width(tilemat)+tiles(2) > 4
    %         set(gcf, 'Color', [1,1,1], 'Position', [100 100 1000 600]);
    %     end

    for i = 1:length(tiles)
        nexttile(tiles(i))
        set(gca,'linewidth',1) % set the axis color
        %         numtiles = length(find(tilemat==tiles(i)));
        %         pbaspect auto
        %         pbaspect(gca,[numtiles 1 1]); % make axes square
        box off
        if length(xticks) > 50
            ticks = xticks;
            set(gca, 'XTick', ticks(1:10:end))
        elseif length(xticks) > 5 % && isnumeric(xticklabels)
            ticks = xticks;
            set(gca, 'XTick', ticks(1:2:end))
        end
        if length(yticks)>50 % && isnumeric(yticklabels)
            ticks = yticks;
            set(gca, 'YTick', ticks(1:10:end))
        elseif length(yticks)>5 % && isnumeric(yticklabels)
            ticks = yticks;
            set(gca, 'YTick', ticks(1:2:end))
        end
        set(gca,'TickDir','out');
        temp = gca;
        temp.Title.FontSize = 10;
        temp.YLabel.FontSize = 10;
        temp.XLabel.FontSize = 10;

    end

else
    set(gca,'linewidth',1) % set the axis color
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
    temp.Title.FontSize = 10;
    temp.YLabel.FontSize = 10;
    temp.XLabel.FontSize = 10;

end


end

