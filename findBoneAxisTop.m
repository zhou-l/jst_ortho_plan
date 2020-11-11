function Iplan = findBoneAxisTop(xRayFname, segMaskBWFile)

I = imread(xRayFname);
dimI = size(I);
if length(dimI)>=3
    I = rgb2gray(I);
end
I = im2double(I);

% load the mask file
M = imread(segMaskBWFile);
% imshow(M);
Mui = im2uint8(M);
Mui = cat(3, Mui, Mui, Mui);
Mui(:,:,1) = 0;
Mui(:,:,2) = 0;



% imshow(I);
% find skeleton
% SK = bwmorph(gpuArray(M), 'skel', Inf);
 SK = bwmorph(M, 'thin', Inf);

% find the contour of the image
% [B, L] = bwboundaries(M);
imshow(labeloverlay(I, SK, 'Transparency', 0));