function img_mat = readFaceImgs(imgDir,imgFmt,outFilePath) 
% Read face images and save them into a cell array
% 
% Inputs:
%   imgDir - The directory storing face images 
%            Provide a relative path if the directory is inside the project 
%            folder.
%            Provide an absolute path if it is outside the project folder.
%   imgFmt - format of the images (e.g., 'jpg', 'tif', etc.)
%   outFilePath - the full path (including file name) of the mat file
%                 storing the cell array

%% Read face images
files = dir(fullfile(imgDir, ['*.' imgFmt]));
% read every face image 
for i = 1:length(files)
    img = imread(fullfile(imgDir,files(i).name));
    img_mat{i,1} = img;
end

%% (Optional) Save the cell array as mat file
if exist('outFilePath', 'var')
    save(outFilePath, 'img_mat')
end

end
