function eBotVert = processBotBoneVertLine(ROIfname, zisRight)
close all;
I = imread(ROIfname);
if nargin < 2
    isRight = true;
else
    isRight = zisRight;
end

dimI = size(I);
if length(dimI)>=3
    I = rgb2gray(I);
end
I = im2double(I);
imshow(I);