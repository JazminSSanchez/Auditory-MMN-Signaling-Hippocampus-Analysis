function [ pathsSpikes tmp] = ReadPCAPaths(folder)
%Reads all the folders saved in the specified path

if nargin < 3
        files = {};
    end
    if nargin < 2
        dirs = {};
    end
    tmp = dir(folder);
    tmp = tmp(~ismember({tmp.name},{'.' '..'}));

    for i = 1:numel(tmp)
        switch tmp(i).isdir
            case 0 % is file
                [pathstr,name,ext] = fileparts(tmp(i).name);
                files = cat(1, files, fullfile(tmp(i).folder, tmp(i).name));
%             case 1 % is dir
%                 dirpath = fullfile(tmp(i).folder, tmp(i).name);
%                 dirs = cat(1, dirs, dirpath);
%                 [dirs, files] = subfiles(dirpath, dirs, files);
        end
        pathsSpikes(i,1) = join([tmp(i).folder,"\",tmp(i).name]);
        pathsSpikes(i,1) = erase(pathsSpikes(i,1) ," ");
    end
end