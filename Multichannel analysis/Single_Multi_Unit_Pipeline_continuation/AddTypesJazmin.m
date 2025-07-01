% Add Oddball stimuli types to bst files
% Just a loop that applies the BST_AddType function to all BST files.
% Adam Hockley 28/3/24
% Changes done by Jazmin Sanchez Sanchez 01/2025

clear
clc

Main_path = 'D:\Jazmin\MultichannelDataTanks\HIP\';
tanks = {'24_023','24_024','24_030'};
Sorter = 'KS4';

for ta = 1:length(tanks)

    % Get electrode positions
    Positions = allfolders([Main_path tanks{ta}]);

    for pos = 1:length(Positions)

        tank_path = [Main_path tanks{ta} '\' Positions{pos}];

        % Use this bit of code to get all Stim blocks in the tank
        Blocks = allfolders(tank_path);
        Blocks = Blocks(contains(Blocks,'STIM'));

        for bl = 1:length(Blocks)

            block_path = [tank_path '\' Blocks{bl}];

            if exist([block_path '\bst_' Sorter '.mat'])

                load([block_path '\bst_' Sorter '.mat']) % load the bst
                bst = BST_AddType(bst);
                save([block_path '\bst_' Sorter '.mat'],'bst') % load the bst

            end
        end
    end
end

%% This function is used to find oddballs & deviants etc in a bst. It will add the info to the bst epoch allowing the trials to be searchable.
% AH 230531
% Only tested on my own stimuli blocks, anyone else planning to use should
% test for correct function.
% Modified by Jazmin Sanchez 2024
% Added recognition of the first 10 tones of each condition, commented the ommission out

function bst = BST_AddType(bst)

% bst.Epocs.Values.typeFirst=[];

if ~ismember('am1_',bst.Epocs.Values.Properties.VariableNames)
    bst.Epocs.Values.am1_(1:height(bst.Epocs.Values)) = 0;
end


allfreqs = unique(bst.Epocs.Values.wfrq);

% Determining frequency pair as two most common frequencies used
for i = 1:length(allfreqs)
    freqscount(i) = sum(bst.Epocs.Values.wfrq == allfreqs(i));
end
freqsloc = find(freqscount == max(freqscount));
freqpair = allfreqs(freqsloc);
if length(freqpair) ~= 2
    error('Did not determine pair of frequencies correctly')
end

allfreqs = [allfreqs allfreqs]; % loop this to help determine cascades
    n = 0;
    m = 0;
    o = 0;
    p = 0;
% Now loop through all bst and assign an epoch of oddball types to each sweep
for i = 1:height(bst.Epocs.Values) % not starting at 0 because silences
    % ODD ASC
    if bst.Epocs.Values.wfrq(i) == 100
        bst.Epocs.Values.type{i} = '';
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100) && (bst.Epocs.Values.wfrq(i+9) == freqpair(2));
        bst.Epocs.Values.type{i} = 'ODDP';
        bst.Epocs.Values.typeFirst{i} = 'ODD-ASC-1-P';
        bst.Epocs.Values.typeFirst{i+1} = 'ODD-ASC-2-P';
        bst.Epocs.Values.typeFirst{i+2} = 'ODD-ASC-3-P';
        bst.Epocs.Values.typeFirst{i+3} = 'ODD-ASC-4-P';
        bst.Epocs.Values.typeFirst{i+4} = 'ODD-ASC-5-P';
        bst.Epocs.Values.typeFirst{i+5} = 'ODD-ASC-6-P';
        bst.Epocs.Values.typeFirst{i+6} = 'ODD-ASC-7-P';
        bst.Epocs.Values.typeFirst{i+7} = 'ODD-ASC-8-P';
        bst.Epocs.Values.typeFirst{i+8} = 'ODD-ASC-9-P';
        bst.Epocs.Values.typeFirst{i+9} = 'ODD-ASC-10-P';
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-1) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-2) == freqpair(1))  && (bst.Epocs.Values.am1_(i) == 0)  && (bst.Epocs.Values.am1_(i-1) == 0);
        bst.Epocs.Values.type{i} = 'ODD-ASC-DEV';
        n = n + 1;
        if n > 40
            bst.Epocs.Values.type{i} = 'ODD-ASC-DEV-P';
        end
        bst.Epocs.Values.type{i-1} = 'ODD-ASC-STD';
        m = m + 1;
        if m > 40
            bst.Epocs.Values.type{i-1} = 'ODD-ASC-STD-P';
        end
    elseif (n < 40) && (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100);
        bst.Epocs.Values.type{i} = 'ODDR';
        bst.Epocs.Values.typeFirst{i} = 'ODD-ASC-1';
        bst.Epocs.Values.typeFirst{i+1} = 'ODD-ASC-2';
        bst.Epocs.Values.typeFirst{i+2} = 'ODD-ASC-3';
        bst.Epocs.Values.typeFirst{i+3} = 'ODD-ASC-4';
        bst.Epocs.Values.typeFirst{i+4} = 'ODD-ASC-5';
        bst.Epocs.Values.typeFirst{i+5} = 'ODD-ASC-6';
        bst.Epocs.Values.typeFirst{i+6} = 'ODD-ASC-7';
        bst.Epocs.Values.typeFirst{i+7} = 'ODD-ASC-8';
        bst.Epocs.Values.typeFirst{i+8} = 'ODD-ASC-9';
        bst.Epocs.Values.typeFirst{i+9} = 'ODD-ASC-10';
        % ODD DES
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100) && (bst.Epocs.Values.wfrq(i+9) == freqpair(1));
        bst.Epocs.Values.type{i} = 'ODDP';
        bst.Epocs.Values.typeFirst{i} = 'ODD-DES-1-P';
        bst.Epocs.Values.typeFirst{i+1} = 'ODD-DES-2-P';
        bst.Epocs.Values.typeFirst{i+2} = 'ODD-DES-3-P';
        bst.Epocs.Values.typeFirst{i+3} = 'ODD-DES-4-P';
        bst.Epocs.Values.typeFirst{i+4} = 'ODD-DES-5-P';
        bst.Epocs.Values.typeFirst{i+5} = 'ODD-DES-6-P';
        bst.Epocs.Values.typeFirst{i+6} = 'ODD-DES-7-P';
        bst.Epocs.Values.typeFirst{i+7} = 'ODD-DES-8-P';
        bst.Epocs.Values.typeFirst{i+8} = 'ODD-DES-9-P';
        bst.Epocs.Values.typeFirst{i+9} = 'ODD-DES-10-P';
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-1) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-2) == freqpair(2))  && (bst.Epocs.Values.am1_(i) == 0)  && (bst.Epocs.Values.am1_(i-1) == 0);
        bst.Epocs.Values.type{i} = 'ODD-DES-DEV';
        o = o + 1;
        if o > 40
            bst.Epocs.Values.type{i} = 'ODD-DES-DEV-P';
        end 
        bst.Epocs.Values.type{i-1} = 'ODD-DES-STD';
        p = p + 1;
        if p > 40
            bst.Epocs.Values.type{i-1} = 'ODD-DES-STD-P';
        end 
    elseif (o < 40) && (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100);
        bst.Epocs.Values.type{i} = 'ODDR';
        bst.Epocs.Values.typeFirst{i} = 'ODD-DES-1';
        bst.Epocs.Values.typeFirst{i+1} = 'ODD-DES-2';
        bst.Epocs.Values.typeFirst{i+2} = 'ODD-DES-3';
        bst.Epocs.Values.typeFirst{i+3} = 'ODD-DES-4';
        bst.Epocs.Values.typeFirst{i+4} = 'ODD-DES-5';
        bst.Epocs.Values.typeFirst{i+5} = 'ODD-DES-6';
        bst.Epocs.Values.typeFirst{i+6} = 'ODD-DES-7';
        bst.Epocs.Values.typeFirst{i+7} = 'ODD-DES-8';
        bst.Epocs.Values.typeFirst{i+8} = 'ODD-DES-9';
        bst.Epocs.Values.typeFirst{i+9} = 'ODD-DES-10';
        % ODD ASC LED
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-1) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-2) == freqpair(1))  && (bst.Epocs.Values.am1_(i) > 0)  && (bst.Epocs.Values.am1_(i-1) == 0);
        bst.Epocs.Values.type{i} = 'ODD-ASC-DEV-LED';
        bst.Epocs.Values.type{i-1} = 'ODD-ASC-STD-LED';
        % ODD DES LED
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-1) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-2) == freqpair(2))  && (bst.Epocs.Values.am1_(i) > 0)  && (bst.Epocs.Values.am1_(i-1) == 0);
        bst.Epocs.Values.type{i} = 'ODD-DES-DEV-LED';
        bst.Epocs.Values.type{i-1} = 'ODD-DES-STD-LED';
        % ODD ASC LED PRE
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-1) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-2) == freqpair(1))  && (bst.Epocs.Values.am1_(i) == 0)  && (bst.Epocs.Values.am1_(i-1) > 0);
        bst.Epocs.Values.type{i} = 'ODD-ASC-DEV-LED-PRE';
        bst.Epocs.Values.type{i-1} = 'ODD-ASC-STD-LED-PRE';
        % ODD DES LED PRE
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (bst.Epocs.Values.wfrq(i-1) == freqpair(2)) && (bst.Epocs.Values.wfrq(i-2) == freqpair(2))  && (bst.Epocs.Values.am1_(i) == 0)  && (bst.Epocs.Values.am1_(i-1) > 0);
        bst.Epocs.Values.type{i} = 'ODD-DES-DEV-LED-PRE';
        bst.Epocs.Values.type{i-1} = 'ODD-DES-STD-LED-PRE';

    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=20) && (i<421)
        bst.Epocs.Values.type{i} = 'ODDR';
    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=460) && (i<861)
        bst.Epocs.Values.type{i} = 'ODDP';
    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=900) && (i<1301)
        bst.Epocs.Values.type{i} = 'ODDR';
    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=1341) && (i<1740)
        bst.Epocs.Values.type{i} = 'ODDP';

    % CASCs
    elseif ( isequal(bst.Epocs.Values.wfrq(i:i+5),allfreqs(freqsloc(1):freqsloc(1)+5)') || isequal(bst.Epocs.Values.wfrq(i-5:i),allfreqs(freqsloc(1)-5:freqsloc(1))') ) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'CASC-ASC-F1';
    elseif ( isequal(bst.Epocs.Values.wfrq(i:i+5),allfreqs(freqsloc(2):freqsloc(2)+5)') || isequal(bst.Epocs.Values.wfrq(i-5:i),allfreqs(freqsloc(2)-5:freqsloc(2))') ) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'CASC-ASC-F2';
    elseif ( isequal(bst.Epocs.Values.wfrq(i:i+5),allfreqs(freqsloc(1):freqsloc(1)+5)') || isequal(bst.Epocs.Values.wfrq(i-5:i),allfreqs(freqsloc(1)-5:freqsloc(1))') ) && (bst.Epocs.Values.am1_(i) > 0);
        bst.Epocs.Values.type{i} = 'CASC-ASC-F1-LED';
    elseif ( isequal(bst.Epocs.Values.wfrq(i:i+5),allfreqs(freqsloc(2):freqsloc(2)+5)') || isequal(bst.Epocs.Values.wfrq(i-5:i),allfreqs(freqsloc(2)-5:freqsloc(2))') ) && (bst.Epocs.Values.am1_(i) > 0);
        bst.Epocs.Values.type{i} = 'CASC-ASC-F2-LED';
    elseif  (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100) &&  ((bst.Epocs.Values.wfrq(i)) < (bst.Epocs.Values.wfrq(i+1))) && ((bst.Epocs.Values.wfrq(i+1)) < (bst.Epocs.Values.wfrq(i+2))) && ((bst.Epocs.Values.wfrq(i+2)) < (bst.Epocs.Values.wfrq(i+3))) && ((bst.Epocs.Values.wfrq(i+3)) < (bst.Epocs.Values.wfrq(i+4)));
        bst.Epocs.Values.type{i} = 'CASCasc';
        bst.Epocs.Values.typeFirst{i} = 'CASC-ASC-1';
        bst.Epocs.Values.typeFirst{i+1} = 'CASC-ASC-2';
        bst.Epocs.Values.typeFirst{i+2} = 'CASC-ASC-3';
        bst.Epocs.Values.typeFirst{i+3} = 'CASC-ASC-4';
        bst.Epocs.Values.typeFirst{i+4} = 'CASC-ASC-5';
        bst.Epocs.Values.typeFirst{i+5} = 'CASC-ASC-6';
        bst.Epocs.Values.typeFirst{i+6} = 'CASC-ASC-7';
        bst.Epocs.Values.typeFirst{i+7} = 'CASC-ASC-8';
        bst.Epocs.Values.typeFirst{i+8} = 'CASC-ASC-9';
        bst.Epocs.Values.typeFirst{i+9} = 'CASC-ASC-10';
    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=2220) && (i<2621)
        bst.Epocs.Values.type{i} = 'CASCasc';

    elseif ( isequal(bst.Epocs.Values.wfrq(i:-1:i-5),allfreqs(freqsloc(1):freqsloc(1)+5)') || isequal(bst.Epocs.Values.wfrq(i+5:-1:i),allfreqs(freqsloc(1)-5:freqsloc(1))') ) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'CASC-DES-F1';
    elseif ( isequal(bst.Epocs.Values.wfrq(i:-1:i-5),allfreqs(freqsloc(2):freqsloc(2)+5)') || isequal(bst.Epocs.Values.wfrq(i+5:-1:i),allfreqs(freqsloc(2)-5:freqsloc(2))') ) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'CASC-DES-F2';
    elseif ( isequal(bst.Epocs.Values.wfrq(i:-1:i-5),allfreqs(freqsloc(1):freqsloc(1)+5)') || isequal(bst.Epocs.Values.wfrq(i+5:-1:i),allfreqs(freqsloc(1)-5:freqsloc(1))') ) && (bst.Epocs.Values.am1_(i) > 0);
        bst.Epocs.Values.type{i} = 'CASC-DES-F1-LED';
    elseif ( isequal(bst.Epocs.Values.wfrq(i:-1:i-5),allfreqs(freqsloc(2):freqsloc(2)+5)') || isequal(bst.Epocs.Values.wfrq(i+5:-1:i),allfreqs(freqsloc(2)-5:freqsloc(2))') ) && (bst.Epocs.Values.am1_(i) > 0);
        bst.Epocs.Values.type{i} = 'CASC-DES-F2-LED';
    elseif  (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100) &&  ((bst.Epocs.Values.wfrq(i)) > (bst.Epocs.Values.wfrq(i+1))) && ((bst.Epocs.Values.wfrq(i+1)) > (bst.Epocs.Values.wfrq(i+2))) && ((bst.Epocs.Values.wfrq(i+2)) > (bst.Epocs.Values.wfrq(i+3))) && ((bst.Epocs.Values.wfrq(i+3)) > (bst.Epocs.Values.wfrq(i+4)));
        bst.Epocs.Values.type{i} = 'CASCdesc';
        bst.Epocs.Values.typeFirst{i} = 'CASC-DES-1';
        bst.Epocs.Values.typeFirst{i+1} = 'CASC-DES-2';
        bst.Epocs.Values.typeFirst{i+2} = 'CASC-DES-3';
        bst.Epocs.Values.typeFirst{i+3} = 'CASC-DES-4';
        bst.Epocs.Values.typeFirst{i+4} = 'CASC-DES-5';
        bst.Epocs.Values.typeFirst{i+5} = 'CASC-DES-6';
        bst.Epocs.Values.typeFirst{i+6} = 'CASC-DES-7';
        bst.Epocs.Values.typeFirst{i+7} = 'CASC-DES-8';
        bst.Epocs.Values.typeFirst{i+8} = 'CASC-DES-9';
        bst.Epocs.Values.typeFirst{i+9} = 'CASC-DES-10';
    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=2660) && (i<3061)
        bst.Epocs.Values.type{i} = 'CASCdesc';

         % MS
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (length(unique(bst.Epocs.Values.wfrq(i-3:i+3))) > 3) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'MS-F1';
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (length(unique(bst.Epocs.Values.wfrq(i-3:i+3))) > 3) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'MS-F2';
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(1)) && (length(unique(bst.Epocs.Values.wfrq(i-3:i+3))) > 3) && (bst.Epocs.Values.am1_(i) > 0);
        bst.Epocs.Values.type{i} = 'MS-F1-LED';
    elseif (bst.Epocs.Values.wfrq(i) == freqpair(2)) && (length(unique(bst.Epocs.Values.wfrq(i-3:i+3))) > 3) && (bst.Epocs.Values.am1_(i) > 0);
        bst.Epocs.Values.type{i} = 'MS-F2-LED';
    elseif (bst.Epocs.Values.wfrq(i-1) == 100) && (bst.Epocs.Values.wfrq(i-2) == 100) && (length(unique(bst.Epocs.Values.wfrq(i-3:i+3))) > 3) && (bst.Epocs.Values.am1_(i) == 0);
        bst.Epocs.Values.type{i} = 'MSall';
        bst.Epocs.Values.typeFirst{i} = 'MS-1';
        bst.Epocs.Values.typeFirst{i+1} = 'MS-2';
        bst.Epocs.Values.typeFirst{i+2} = 'MS-3';
        bst.Epocs.Values.typeFirst{i+3} = 'MS-4';
        bst.Epocs.Values.typeFirst{i+4} = 'MS-5';
        bst.Epocs.Values.typeFirst{i+5} = 'MS-6';
        bst.Epocs.Values.typeFirst{i+6} = 'MS-7';
        bst.Epocs.Values.typeFirst{i+7} = 'MS-8';
        bst.Epocs.Values.typeFirst{i+8} = 'MS-9';
        bst.Epocs.Values.typeFirst{i+9} = 'MS-10';
    elseif (bst.Epocs.Values.wfrq(i) ~= 100) &&  (i>=1780) && (i<2181)
        bst.Epocs.Values.type{i} = 'MSall';
        % %% OMISSIONS
        % for i = 2:height(data.epocs.Wfrq.data)-1 % not starting at 0 because silences
        % 
        %     if (data.Epocs.Values.Us1_(i) == 0) && (data.Epocs.Values.Us1_(i-1) > 0) && (data.Epocs.Values.Us1_(i+1) > 0) && (data.Epocs.Values.wfrq(i-1) == freqpair(1)) && (data.Epocs.Values.am1_(i) == 0)
        %         output{i} = 'OM-F1';
        %     elseif (data.Epocs.Values.Us1_(i) == 0) && (data.Epocs.Values.Us1_(i-1) > 0) && (data.Epocs.Values.Us1_(i+1) > 0) && (data.Epocs.Values.wfrq(i-1) == freqpair(2)) && (data.Epocs.Values.am1_(i) == 0)
        %         output{i} = 'OM-F2';
        %     % elseif (data.Epocs.Values.Us1_(i) == 0) && (data.Epocs.Values.Us1_(i-1) > 0) && (data.Epocs.Values.Us1_(i+1) > 0) && (data.Epocs.Values.wfrq(i-1) == freqpair(1)) && (data.Epocs.Values.am1_(i) > 0)
        %     %     output{i} = 'OM-F1-LED';
        %     % elseif (data.Epocs.Values.Us1_(i) == 0) && (data.Epocs.Values.Us1_(i-1) > 0) && (data.Epocs.Values.Us1_(i+1) > 0) && (data.Epocs.Values.wfrq(i-1) == freqpair(2)) && (data.Epocs.Values.am1_(i) > 0)
        %     %     output{i} = 'OM-F2-LED';
        %     end
        % end
        % %%
        % % STDs
        % for i = 2:height(data.epocs.Wfrq.data)-1 % not starting at 0 because silences
        % 
        %     if (data.Epocs.Values.wfrq(i) == freqpair(1)) &&  (data.Epocs.Values.wfrq(i-1) == freqpair(1)) &&  (data.Epocs.Values.wfrq(i+1) == freqpair(1)) && (data.Epocs.Values.am1_(i) > 0)
        %         output{i} = 'STDS-F1-LED';
        %     elseif (data.Epocs.Values.wfrq(i) == freqpair(2)) &&  (data.Epocs.Values.wfrq(i-1) == freqpair(2)) &&  (data.Epocs.Values.wfrq(i+1) == freqpair(2)) && (data.Epocs.Values.am1_(i) > 0)
        %         output{i} = 'STDS-F2-LED';
        %     end

            %     if (data.Epocs.Values.Us1_(i) == 0) && (data.Epocs.Values.Us1_(i-1) > 0) && (data.Epocs.Values.Us1_(i+1) > 0) && (data.Epocs.Values.wfrq(i-1) == freqpair(1)) && (data.Epocs.Values.am1_(i) > 0)
            %         output{i} = 'STDS-F1-LED';
            %     elseif (data.Epocs.Values.Us1_(i) == 0) && (data.Epocs.Values.Us1_(i-1) > 0) && (data.Epocs.Values.Us1_(i+1) > 0) && (data.Epocs.Values.wfrq(i-1) == freqpair(2)) && (data.Epocs.Values.am1_(i) > 0)
            %         output{i} = 'STDS-F2-LED';
            %     end
        end
    end
end

% Here test if even nnumber of each condition was typed
types = {'ODD-ASC-DEV','ODD-ASC-STD','ODD-DES-DEV','ODD-DES-STD','ODD-ASC-DEV-P',...
    'ODD-ASC-STD-P','ODD-DES-DEV-P','ODD-DES-STD-P','CASC-ASC-F1',...
    'CASC-ASC-F2','CASC-DES-F1','CASC-DES-F2','MS-F1','MS-F2'};

clear trials

for t = 1:length(types)
    SelectedTrials = BST_TS3(bst,'type',types{t}); % get selected trials of freqs, levels, estim etc... any epoch
    trials(t) = length(SelectedTrials);
end

if length(unique(trials)) > 1
    warning([bst.tank(end-6:end) '\' bst.Block ': Uneven number of typed conditions'])
end
% end

%%
function folders = allfolders(directory)

folders = dir(directory);
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
folders = folders(dirFlags);
folders = {folders.name};

end
