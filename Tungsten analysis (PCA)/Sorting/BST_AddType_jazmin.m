function bst = BST_AddType_jazmin(bst)

%% This function is used to find oddballs & deviants etc in a bst. It will add the info to the bst epoch allowing the trials to be searchable.
% AH 230531
% Only tested on my own stimuli blocks, anyone else planning to use should
% test for correct function.
% Modified by Jazmin Sanchez 231201
for k = 1:height(bst.Epocs.Values.freq)
    s = 0;
    s = sprintf('%.0f', bst.Epocs.Values.freq(k));
    bst.Epocs.Values.freq(k) = str2double( s );
end
allfreqs = unique(bst.Epocs.Values.freq);

% Determining frequency pair as two most common frequencies used
for i = 1:length(allfreqs)
    freqscount(i) = sum(bst.Epocs.Values.freq == allfreqs(i));
end
freqsloc = find(freqscount == max(freqscount));
freqpair = allfreqs(freqsloc);
if length(freqpair) ~= 2
    error('Did not determine pair of frequencies correctly')
end

allfreqs = [allfreqs allfreqs]; % loop this to help determine cascades

% Now loop through all bst and assign an epoch of oddball types to each sweep
for i = 1:height(bst.Epocs.Values) % not starting at 0 because silences
    % ODD ASC
    if bst.Epocs.Values.freq(i) == 100
        bst.Epocs.Values.type{i} = '';
    elseif (contains(bst.Epocs.Values.Block(i),'ODD') == 1) && (bst.Epocs.Values.freq(i) == freqpair(2)) && (bst.Epocs.Values.freq(i-1) == freqpair(1)) && (bst.Epocs.Values.freq(i-2) == freqpair(1));
        bst.Epocs.Values.type{i} = 'ODD-ASC-DEV';
        bst.Epocs.Values.type{i-1} = 'ODD-ASC-STD';
        % ODD DES
    elseif (contains(bst.Epocs.Values.Block(i),'ODD') == 1) && (bst.Epocs.Values.freq(i) == freqpair(1)) && (bst.Epocs.Values.freq(i-1) == freqpair(2)) && (bst.Epocs.Values.freq(i-2) == freqpair(2));
        bst.Epocs.Values.type{i} = 'ODD-DES-DEV';
        bst.Epocs.Values.type{i-1} = 'ODD-DES-STD';

    % CASCs

    elseif (contains(bst.Epocs.Values.Block(i),'CASCASC') == 1) && (bst.Epocs.Values.freq(i) == freqpair(1)) && (bst.Epocs.Values.freq(i+1) == freqpair(2));
        bst.Epocs.Values.type{i} = 'CASC-ASC-F1';
    elseif (contains(bst.Epocs.Values.Block(i),'CASCASC') == 1) && (bst.Epocs.Values.freq(i) == freqpair(2)) && (bst.Epocs.Values.freq(i-1) == freqpair(1));
        bst.Epocs.Values.type{i} = 'CASC-ASC-F2';


    elseif (contains(bst.Epocs.Values.Block(i),'CASCDESC') == 1) && (bst.Epocs.Values.freq(i) == freqpair(1)) && (bst.Epocs.Values.freq(i-1) == freqpair(2));
        bst.Epocs.Values.type{i} = 'CASC-DES-F1';
    elseif (contains(bst.Epocs.Values.Block(i),'CASCDESC') == 1) && (bst.Epocs.Values.freq(i) == freqpair(2)) && (bst.Epocs.Values.freq(i+1) == freqpair(1));
        bst.Epocs.Values.type{i} = 'CASC-DES-F2';


    end
end
bst.Epocs.TSOn.type = bst.Epocs.TSOn.freq;
bst.Epocs.TSOff.type = bst.Epocs.TSOff.freq;
bst.NTrials = length(bst.Epocs.Values.freq);
