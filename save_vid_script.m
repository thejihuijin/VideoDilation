clear; clc;
basePath = 'data/Reddit_Videos/';
outpath = strcat(basePath,'Dilated/');
myFiles = dir(fullfile(basePath,'*.mp4'));
savevidtime = tic;
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  if ~strcmp(baseFileName(1:3),'rs_')
    filename = fullfile(basePath, baseFileName);
    fprintf(1, 'Now reading %s\n', filename);
    outputfilename = fullfile(outpath, strcat('dilated_',baseFileName));
    dilate_video;
  end
  
end
toc(savevidtime);