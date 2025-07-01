function [Spikes] = BST_GS(bst,trials,rastertype)

Spikes=bst.Spikes.(['Raster' rastertype])(ismember(bst.Spikes.TrialIdx,trials));

% NOT NEEDED BUT KEPT JUST IN CASE (old junk XM code from sst days)
% [~,sep] = unique(bst.Spikes.BlockIdx);
% sep=[sep;length(bst.Spikes.BlockIdx)+1];
% for i=2:length(sep)-1
%     tank_n = bst.Spikes.BlockIdx(sep(i));
%     trial_prev = max(find(bst.Epocs.Values.tind==(tank_n-1)));
% 
%     bst.Spikes.TrialIdx(sep(i):sep(i+1)-1) = bst.Spikes.TrialIdx(sep(i):sep(i+1)-1) + trial_prev;
% end

