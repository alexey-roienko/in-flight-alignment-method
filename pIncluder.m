
cd(['logs' '\' logFolderName '\']);
logsFolder = pwd;

cd(['..' '\' '..' '\' 'classes' '\']);
classesFolder = pwd;

cd(['..' '\' 'plugins' '\']);
pluginsFolder = pwd;

cd(['..' '\' 'subroutings' '\']);
subrFolder = pwd;

cd(curDir);

addpath(classesFolder, logsFolder, pluginsFolder, subrFolder);
clear classesFolder pluginsFolder subrFolder