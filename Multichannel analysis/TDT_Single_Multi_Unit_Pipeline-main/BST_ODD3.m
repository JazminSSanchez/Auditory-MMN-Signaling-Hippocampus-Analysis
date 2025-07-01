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

