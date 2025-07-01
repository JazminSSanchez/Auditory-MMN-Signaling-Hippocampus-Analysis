function sb_cell = build_sb_jazmin(path, save_dir, rfs_user, varargin)

% path={'D:\Tanks\AdamTest-230118-173319'}
% rfs_user={'fullFRA'};
% save_dir=''D:\Tanks\tempSB' %temporary superblock when just looking at RFs during experiment
% varargin, put anything to disable remove electrial
%
% This code allow users to construct SUPERBLOCK (sb_cell: cell array with
% only one field), which is the first step for spika analysis. superblocks
% cell is then fed into bbst in the next step.
%   
% Some of the stuff in here is obscelete from previous labs etc. 
%
% - AH
%

keep_artifact=true;
if ~isempty(varargin)
    keep_artifact=true;
end

% check input consistency
if (isstr(path)&size(rfs_user,1)~=1)|(iscell(path)&size(rfs_user,1)==1)
    error('check tank path and block list')
end
if size(rfs_user,1)==2&max(rfs_user(2,:))>length(cellstr(path))
    error('block index inconsistent')
end

save_dir=fullfile(save_dir,'superblock');
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end


if isstr(path)
    sb_out = build_rfblock_named(path, rfs_user, keep_artifact);
elseif iscell(path)
    for i=1:length(path)
        path_str=path{i};
        block=rfs_user(1,rfs_user(2,:)==i);
        try
            sb = build_rfblock_simple(path_str, block,keep_artifact);
        catch
            error([path_str ' does not exist']);
        end
        if i==1
            sb_out=[array2table(ones(height(sb),1)*i,'variablenames',{'tank'}) sb];
        else
            sb_out=[sb_out;array2table(ones(height(sb),1)*i,'variablenames',{'tank'}) sb];
        end
        fprintf('Tank-%d (of %d) compiled.\n',i,length(path));
        %             fprintf(fid,'%s: Tank-%d (of %d) compiled.\r\n',datestr(now),i,length(path));
    end
else
    error('Invalid tank path format');
end
sb_cell={};
sb_cell{1}=sb_out;
save_separate_channels(save_dir,sb_out);

%     fclose(fid);
end


function save_separate_channels(save_dir, sb_out)

fprintf('Saving...\n')
% fprintf(fid,'%s: Saving...\r\n',datestr(now));

max_ch=max(sb_out.chan);
for channel=1:max_ch
    sb = sb_out(sb_out.chan==channel,:);
    fname=fullfile(save_dir,['channel_' num2str(channel) '.mat']);
    save(fname,'sb','-v7.3')
end
fprintf('n = %d channels saved.\n',max_ch)
% 	fprintf(fid,'%s: n = %d channels saved.\r\n',datestr(now),max_ch);

end



function superblocks = build_rfblock_named(path, rfs_user, keep_artifact)

% see BUILD_RFBLOCK & build_rfblock_simple
% This code is altered to allow naming of blocks, not simply block-1 etc. 
% 
% - AH


for i_block=1:length(rfs_user)
    
    blockStr = rfs_user{i_block};
        blockN=i_block;

    % Loading in the data, if theres already a mat file (from
    % rethresholding) then use that, else do TDT2mat
    if exist([path '\' blockStr '_thresh.mat'])
        load([path '\' blockStr '_thresh.mat'])
        fprintf('Using rethresholded data\n');
    else
        try
%             data=TDT2mat(path,blockStr,'Type',[2 3],'Verbose',0);
            data=TDTbin2mat([path '/' blockStr],'Type',[2 3],'Verbose',0);
            
        catch
            error([path '\' blockStr ' does not exist']);
        end
    end
    
    if ~isempty(data.snips)
        
        spike_var = fieldnames(data.snips);
        if strcmp(spike_var{1},'CSPK')
            fprintf('Legacy Tank (OpenEx)\n');
        end
        block=ones(length(data.snips.(spike_var{1}).chan),1)*blockN;
        chan=data.snips.(spike_var{1}).chan;
        ts=data.snips.(spike_var{1}).ts;
        waves=double(data.snips.(spike_var{1}).data);
        sortc=zeros(length(data.snips.(spike_var{1}).chan),1);
        clear BlockIdx
        BlockIdx(1:length(data.snips.(spike_var{1}).chan),1) = i_block;
        
        partList=unique(data.epocs.Freq.data);
        part=zeros(length(data.snips.(spike_var{1}).chan),1);
        for i_p=1:length(partList)
            idx=find(data.epocs.Freq.data==partList(i_p));
            t_start=data.epocs.Freq.onset(idx(1));
            t_end=data.epocs.Freq.offset(idx(end));
            ts_idx=find(data.snips.(spike_var{1}).ts>=t_start&data.snips.(spike_var{1}).ts<=t_end);
            part(ts_idx)=partList(i_p);
        end
        
        SB_com=table(block,BlockIdx,chan,ts,sortc,waves,part);
        SB_com(SB_com.part==0,:)=[];
        
        epocs = data.epocs;
        if ~keep_artifact
            e_idx = remove_artifact(SB_com.ts,epocs);
            SB_com(logical(e_idx),:)=[];
        end
    else
        SB_com=table([],[],[],[],[],[],'variablenames',...
            {'block','chan','ts','sortc','waves','part'});
    end
    fprintf('%s compiled. (%d/%d)\n',blockStr,i_block,length(rfs_user))
    %     fprintf(fid,'%s: %s compiled. (%d/%d)\r\n',datestr(now),blockStr,i_block,length(rfs_user));
    
    if i_block==1
        superblocks=SB_com;
    else
        superblocks=[superblocks;SB_com];
    end
end


end

function superblocks = build_rfblock_simple(path, rfs_user, keep_artifact)

% see BUILD_RFBLOCK
% This code modified uneccessary considerations for identifying RFs and
% whatnot. This newer version only requests user for a list of blocks in
% array format and builds a superblock, regardless of which block has
% RF(s).
%
% - CW


for i_block=1:length(rfs_user)
    
    blockN=rfs_user(i_block);
    blockStr=['Block-' num2str(blockN)];
    
    % Loading in the data, if theres already a mat file (from
    % rethresholding) then use that, else do TDT2mat
    if exist([path '\' blockStr '_thresh.mat'])
        load([path '\' blockStr '_thresh.mat'])
        fprintf('Using rethresholded data\n');
    else
        try
%             data=TDT2mat(path,blockStr,'Type',[2 3],'Verbose',0);
            data=TDTbin2mat([path '/' blockStr],'Type',[2 3],'Verbose',0);
            
        catch
            error([path '\' blockStr ' does not exist']);
        end
    end
    
    if ~isempty(data.snips)
        
        spike_var = fieldnames(data.snips);
        if strcmp(spike_var{1},'CSPK')
            fprintf('Legacy Tank (OpenEx)\n');
        end
        block=ones(length(data.snips.(spike_var{1}).chan),1)*blockN;
        chan=data.snips.(spike_var{1}).chan;
        ts=data.snips.(spike_var{1}).ts;
        waves=double(data.snips.(spike_var{1}).data);
        sortc=zeros(length(data.snips.(spike_var{1}).chan),1);
        clear BlockIdx
        BlockIdx(1:length(data.snips.(spike_var{1}).chan),1) = i_block;
        
        partList=unique(data.epocs.FInd.data);
        part=zeros(length(data.snips.(spike_var{1}).chan),1);
        for i_p=1:length(partList)
            idx=find(data.epocs.FInd.data==partList(i_p));
            t_start=data.epocs.FInd.onset(idx(1));
            t_end=data.epocs.FInd.offset(idx(end));
            ts_idx=find(data.snips.(spike_var{1}).ts>=t_start&data.snips.(spike_var{1}).ts<=t_end);
            part(ts_idx)=partList(i_p);
        end
        
        SB_com=table(block,BlockIdx,chan,ts,sortc,waves,part);
        SB_com(SB_com.part==0,:)=[];
        
        epocs = data.epocs;
        if ~keep_artifact
            e_idx = remove_artifact(SB_com.ts,epocs);
            SB_com(logical(e_idx),:)=[];
        end
    else
        SB_com=table([],[],[],[],[],[],'variablenames',...
            {'block','chan','ts','sortc','waves','part'});
    end
    fprintf('%s compiled. (%d/%d)\n',blockStr,i_block,length(rfs_user))
    %     fprintf(fid,'%s: %s compiled. (%d/%d)\r\n',datestr(now),blockStr,i_block,length(rfs_user));
    
    if i_block==1
        superblocks=SB_com;
    else
        superblocks=[superblocks;SB_com];
    end
end


end

function idx_tot = remove_artifact(ts,epocs)

idx_tot = zeros(length(ts),1);
if ~isfield(epocs,'EAmp')|~isfield(epocs,'ETyp')
    return;
end

idx_e = find(epocs.EAmp.data~=0);
if isempty(idx_e)
    return;
end

for i=1:length(idx_e)
    e_on=epocs.ETyp.onset(idx_e(i));
    e_off=epocs.ETyp.offset(idx_e(i));
    idx_tot = idx_tot + double(ts>=e_on&ts<=e_off);
end

fprintf('Removed artifact.\n')
% fprintf(fid,'%s: Removed artifact.\r\n',datestr(now));


end