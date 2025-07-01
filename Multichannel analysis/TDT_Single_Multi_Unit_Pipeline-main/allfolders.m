function folders = allfolders(directory)

folders = dir(directory);
dirFlags = [folders.isdir] & ~strcmp({folders.name},'.') & ~strcmp({folders.name},'..');
folders = folders(dirFlags);
folders = {folders.name};

end