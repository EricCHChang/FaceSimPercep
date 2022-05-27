% Store face images into a mat file (ims_new.mat)

%% 
rootDir = cd;
imgDir = fullfile(rootDir,'stims');
outFilePath = fullfile(rootDir, 'ims_new.mat');
imgFmt = 'tif';

img_mat = readFaceImgs(imgDir,imgFmt,outFilePath); 