function eTopVert = processTopBoneVertLine(ROIfname, zisRight)
% close all;
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
[gradThresh,numIter] = imdiffuseest(I,'ConductionMethod','quadratic');
tmp = imdiffusefilt(I,'ConductionMethod','quadratic', ...
    'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
tmp = impyramid(tmp,'reduce');
tmp = imdiffusefilt(tmp);
[Gmag, Gdir]=imgradient(tmp);
[BWE,edgeThres] = edge(tmp,'Canny');
[BWE,edgeThres2] = edge(tmp,'Canny', [1.5 *edgeThres(2), 2 * edgeThres(2)]);
figure, 
imshowpair(tmp, Gmag, 'montage');
figure;
imshowpair(tmp, BWE, 'montage');
%% extract high contrast boundary of the top part
% normalize values
Gmag = rescale(Gmag, 0, 1);
tmp = rescale(tmp, 0, 1);

gradThresUpperBone = 0.25;
valThresUpperBone = 0.8;

uidx1 = Gmag(:,:) > gradThresUpperBone; % strong boundaries
uidx2 = tmp(:,:) > valThresUpperBone; 

idx = uidx1; %uidx2 & uidx1 ;

tmpX = zeros(size(tmp));
tmpX(idx) = tmp(idx);

figure, imshowpair(tmp, tmpX, 'montage');
StrongE = imbinarize(tmpX);
BWE = BWE & StrongE; 

BWE = bwmorph(BWE,'bridge', inf);
BWE = bwmorph(BWE, 'spur');
imshowpair(StrongE, BWE, 'montage');
%% the bone edges are satisfactory
CC = bwconncomp(BWE) 
% extract the inner bone
% get long-enough CCs (4 longest CCs)
numPixels = cellfun(@numel,CC.PixelIdxList);
[nB,nI] = sort(numPixels,'descend');
minSy = realmax;
% get the top most CC    
inInBoneCC = 0; 
inOutBoneCC = 0;
topCCNum = min(4,numel(CC.PixelIdxList));

% keep only the top 4 CCs
skel = false(size(tmp));
muXlist = [];%zeros(length(CC.PixelIdxList),1)
for i=1:topCCNum
    pxList = CC.PixelIdxList{nI(i)};
    [ptSy,ptSx] = ind2sub(size(tmp),pxList);
    mSx = mean(ptSx);
    muXlist(end+1,:) = [mSx, nI(i)];
%     skel(CC.PixelIdxList{i}) = true;
end
[muSB,muSI] = sort(muXlist,'ascend');
    
if isRight
   inInBoneCC = muSB(1,2);%left most long CC
   inOutBoneCC = muSB(2,2); %second left most long CC
else
   inInBoneCC = muSB(4,2);%right most long CC
   inOutBoneCC = muSB(3,2); %second right most long CC
end

skel(CC.PixelIdxList{inInBoneCC}) = true;
skel(CC.PixelIdxList{inOutBoneCC}) = true;
figure, imshow(skel);
%% compute the vertical line
pxInInBCC = CC.PixelIdxList{inInBoneCC};
[pxInIny, pxInInx] = ind2sub(size(tmp), pxInInBCC);
pxInIn = [pxInIny, pxInInx];
% sort by y
pxInIn = sort(pxInIn, 'ascend');
pxInOutBCC = CC.PixelIdxList{inOutBoneCC};
[pxInOuty, pxInOutx] = ind2sub(size(tmp), pxInOutBCC);
pxInOut = [pxInOuty, pxInOutx];
pxInOut = sort(pxInOut, 'ascend');
% scale back to the original image
ptS = [2 * size(tmp,1),0.5 * 2  * (pxInIn(1,2) + pxInOut(1,2))];
ptT = [1,  0.5 * 2 * (pxInIn(end,2)+pxInOut(end,2))];

figure, imshow(I); hold on;
eTopVert = [ptS, ptT];
plot([eTopVert(2),eTopVert(4)], [eTopVert(1),eTopVert(3)],'r');
hold off;
