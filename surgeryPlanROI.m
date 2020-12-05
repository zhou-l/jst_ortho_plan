% stretchAlpha is the angle (in degrees) of stretching
function Iplan = surgeryPlanROI( xRayInputfName, upperBoneROIfName, lowerBoneROIfName, stretchAlpha, isRightInImg)
if nargin < 5
    zStretchAlpha = 10; % stretch 10 degrees
    zIsRightInImg = true;
else
    zStretchAlpha = stretchAlpha;
    if nargin < 7
        zIsRightInImg = true;
    else
        zIsRightInImg = isRightInImg;
    end
end
close all;
figure;
Ifull = imread(xRayInputfName);
dimI = size(Ifull);
if length(dimI)>=3
    Ifull = rgb2gray(Ifull);
end
Ifull = im2double(Ifull);

hI = imshow(Ifull);hold on;
I_xRayInput = Ifull;
I_upperBoneROI = imread(upperBoneROIfName);
I_lowerBoneROI = imread(lowerBoneROIfName);

[Gmag, Gdir] = imgradient(Ifull,'sobel');
imshow(Gmag);
% alpha(hI, 0.7);
hold on;
isLower = false;
[hAxUpper, vAxUpper, boundUpper, kUpper] = findUpperBoneAxis(I_xRayInput, I_upperBoneROI);
hold on;
isLower = true;
[hAxLower, vAxLower, boundLower, kLower] = findLowerBoneAxis(I_xRayInput, I_lowerBoneROI);

scale = 600;