function eTop = preprocessROItop(ROIfname, zisRight)
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
% tmp = imdiffusefilt(I);

%use scale space???

[gradThresh,numIter] = imdiffuseest(I,'ConductionMethod','quadratic');
tmp = imdiffusefilt(I,'ConductionMethod','quadratic', ...
    'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
tmp = impyramid(tmp,'reduce');

tmp = imdiffusefilt(tmp);
tmp = imdiffusefilt(tmp);
[Gmag, Gdir]=imgradient(tmp);
imshowpair(tmp, Gmag, 'montage');

BWE = edge(tmp, 'Canny');
imshowpair(tmp, BWE, 'montage');
CC = bwconncomp(BWE);
% get CCs that are positioned on top
minLen = 15;
muYlist = [];%zeros(length(CC.PixelIdxList),1);
for i = 1:length(CC.PixelIdxList)
    pxList = CC.PixelIdxList{i};
    if numel(pxList) > minLen
        [ptSy,ptSx] = ind2sub(size(tmp),pxList);
        muY = mean(ptSy);
        muYlist(end+1,:) = muY;
    end
end

[muSB,muSI] = sort(muYlist);
bwUpper = false(size(BWE));
for i = 1:10
    pxList = CC.PixelIdxList{muSI(i)};
    [ptSy,ptSx] = ind2sub(size(tmp),pxList);
    bwUpper(pxList) = true;
end
imshowpair(BWE, bwUpper, 'montage');
% normalize values
Gmag = rescale(Gmag, 0, 1);
tmp = rescale(tmp, 0, 1);


% extract regions with Low Gmag and High intensity
Iupper = zeros(size(tmp));
Ilower = zeros(size(tmp));

gradThresUpperBone = 0.2;
valThresUpperBone = 0.7;

uidx1 = Gmag(:,:) < gradThresUpperBone;
uidx2 = tmp(:,:) > valThresUpperBone; 

idx = uidx2 & uidx1 ;
% idx = uidx2;
bwUpper(idx) = BWE(idx);

% show the thresholded image and ask the user to draw a point
imshowpair(tmp, bwUpper, 'montage');

bwUpper = bwmorph(bwUpper, 'bridge'); % bridge small gaps in the BW
bwUpper = bwmorph(bwUpper, 'spur');

CC = bwconncomp(bwUpper) 
figure;
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
% bwUpper(CC.PixelIdxList{idx}) = 0;
imshow(bwUpper);
% h = drawpoint;
[xi, yi, but] = ginput(1);
xi = int16(xi);
yi = int16(yi);
ptIdx = int16((xi-1)*size(tmp,1)+yi);
ptIdx2 = sub2ind(size(tmp), yi, xi);
mainCC = 0;
for i = 1:length(CC.PixelIdxList)
    pxList = CC.PixelIdxList{i};
    if isempty(find(pxList == ptIdx, 1))
        bwUpper(CC.PixelIdxList{i})= 0;
        continue;
    else
        disp('CC found!');
        mainCC = i;
    end
end
if mainCC < 1
    warning('top boundary not found!');
    return;
end
% trace end points without cutBranch
lastVecL = [0 1];
simpK = [];
% get end points
extendTopPt = [realmax,-realmax];
extendEndPt = [-realmax,-realmax];
%Bd = ind2sub(size(Bskel),skelInd);
skelInd = CC.PixelIdxList{mainCC};
for i = 1:numel(CC.PixelIdxList{mainCC})-1
    j = skelInd(i); j2 = skelInd(i+1);
    [endPtY, endPtX] = ind2sub(size(tmp), j);
    [sPtY,sPtX] = ind2sub(size(tmp), j2);
    line = [sPtY,sPtX,endPtY,endPtX];
    vecL = [line(3)-line(1),line(4)-line(2)];
    % get top and end pts
    if endPtY < extendTopPt(1)  % y is flipped in Matlab-- top is min
         extendTopPt(1) = endPtY;
         extendTopPt(2) = endPtX;
    end
    
    if endPtX > extendEndPt(2)
        extendEndPt(2) = endPtX;
        extendEndPt(1) = endPtY;
    end
    
    if norm(vecL) > 1e-4 
        vecL = vecL / norm(vecL);
    else
        vecL = [0 0];
    end
    % angle difference is small
    if dot(lastVecL,vecL)>cos(pi/30)
        if i == numel(skelInd)-1
            simpK(end+1,:) = j2;
        end
        continue;     
    else
        simpK(end+1,:) = j;
        if i == numel(skelInd)-1
            simpK(end+1,:) = j2;
        end
        lastVecL = vecL;     
    end
end

figure; imshow(tmp); hold on;
[simpPtY,simpPtX] = ind2sub(size(tmp), simpK);
plot(simpPtX,simpPtY, 'LineWidth', 3); 

% scale back to the original image
extendTopPt = extendTopPt .* 2;
extendEndPt = extendEndPt .* 2;
figure; imshow(I); hold on;
eTop = [extendTopPt, extendEndPt];
plot([extendTopPt(2),extendEndPt(2)], [extendTopPt(1),extendEndPt(1)], 'r');
hold off;



