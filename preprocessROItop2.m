function eTop = preprocessROItop2(ROIfname, zisRight)
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

[Gmag, Gdir]=imgradient(tmp);

% normalize values
Gmag = rescale(Gmag, 0, 1);
tmp = rescale(tmp, 0, 1);

%% extract high intensity boundary of the top part
gradThresUpperBone = 0.2;
valThresUpperBone = 0.75;

uidx1 = Gmag(:,:) < gradThresUpperBone;
uidx2 = tmp(:,:) > valThresUpperBone; 

% idx = uidx2 & uidx1 ;
idx = uidx2; 

tmpX = zeros(size(tmp));
tmpX(idx) = tmp(idx);
[BWEo,threshOut]  = edge(tmp, 'Canny');
[BWEo, thresOut2] = edge(tmp, 'Canny', [threshOut(2)*0.8, threshOut(2) * 1.5]);
figure, imshow(BWEo);


imshowpair(tmp, Gmag, 'montage');
figure, imshowpair(tmpX, tmp, 'montage');
[BWE,threshOut]  = edge(tmpX, 'Canny');
figure, imshow(BWE);
% bwUpper = bwmorph(BWE, 'bridge',5); % bridge small gaps in the BW
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
%% fetch top 10 mean y CCs
topCCNum = min([10,length(muYlist),numel(CC.PixelIdxList)]);
for i = 1:topCCNum
    pxList = CC.PixelIdxList{muSI(i)};
    [ptSy,ptSx] = ind2sub(size(tmp),pxList);
    bwUpper(pxList) = true;
end
imshowpair(BWE, bwUpper, 'montage');
% idx = uidx2;
% bwUpper(idx) = BWE(idx);

% show the thresholded image and ask the user to draw a point
imshowpair(tmp, bwUpper, 'montage');

% bwUpper = bwmorph(bwUpper, 'bridge'); % bridge small gaps in the BW
% bwUpper = bwmorph(bwUpper, 'fill', Inf); 
% bwUpper = bwmorph(bwUpper, 'spur');
% bwUpper = bwmorph(bwUpper, 'skel', Inf);
CC = bwconncomp(bwUpper) 
figure; imshow(bwUpper);
numPixels = cellfun(@numel,CC.PixelIdxList);
[nB,nI] = sort(numPixels,'descend');
minSy = realmax;
%% get the top most CC
for i = 1:topCCNum
    pxList = CC.PixelIdxList{nI(i)};
   [ptSy,ptSx] = ind2sub(size(tmp),pxList);
    mSy = mean(ptSy);
   if mSy < minSy
       minSy = mSy;
       mainCC = nI(i);
       endPt = [ptSy(end),ptSx(end)];
   end
end
if mainCC < 1
    warning('top boundary not found!');
    return;
end
skel = false(size(tmp));
skel(CC.PixelIdxList{mainCC}) = true;
figure, imshow(skel); 
%% and search for the vertical CC
bwV = edge(tmp, 'Canny');
VCC = bwconncomp(bwV) ;
vertCC = 0;
minLenThres = 20;

minXDist = realmax;
mainPxList = CC.PixelIdxList{mainCC};
   [ptSy,ptSx] = ind2sub(size(tmp),pxList);
    mainPxListMeanY = mean(ptSy);
for i = 1:length(VCC.PixelIdxList)
    pxList = VCC.PixelIdxList{i};
    if numel(pxList) <= minLenThres
        continue;
    end
    
    [ptSy,ptSx] = ind2sub(size(tmp),pxList(1));
    [ptTy,ptTx] = ind2sub(size(tmp),pxList(end));
    % sort Y from small to large
    if ptSy > ptTy
        ys = ptSy;
        ptSy = ptTy;
        ptTy = ys;
        
        xs = ptSx;
        ptSx = ptTx;
        ptTx = xs;
        
    end
    dir = [ptTx - ptSx, ptTy - ptSy];
    dir = normalize(dir,'norm',2);
    %% NOTE: If it's right leg

    if isRight
        theta = acos(dot(dir,[-1,0]));
    else
        theta = acos(dot(dir,[1,0]));
    end
    if theta - pi/2 < 0 && -theta + pi/2 < pi/4 % The vertical edge should go from bottom right to top left
%         skel(VCC.PixelIdxList{i}) = true;
        % and the vertical range of the vertical edge should cover the main CC
        if mainPxListMeanY >= ptSy && mainPxListMeanY <= ptTy
            dist =  abs(min(norm(endPt - [ptSy,ptSx]), norm(endPt-[ptTy,ptTx])));
            if dist < minXDist
                minXDist = dist;
                vertCC = i;

            end
        end
    end 
end
skel(VCC.PixelIdxList{vertCC}) = true;
figure, imshow(skel); hold on;
lastVecL = [0 1];
simpK = [];
% get end points of the horizontal boundary
extendTopPt = [realmax,-realmax];  % top left
extendEndPt = [-realmax,-realmax]; % bottom right
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
    
    if endPtX > extendEndPt(2) && endPtX > extendTopPt(2) && endPtY > extendEndPt(1)
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


lastVecL = [0 1];
simpK = [];
% get end points
extendBotPtVert = [-realmax,-realmax];
extendEndPtVert = [realmax,-realmax];
%Bd = ind2sub(size(Bskel),skelInd);
skelInd = VCC.PixelIdxList{vertCC};
for i = 1:numel(skelInd)-1
    j = skelInd(i); j2 = skelInd(i+1);
    [endPtY, endPtX] = ind2sub(size(tmp), j);
    [sPtY,sPtX] = ind2sub(size(tmp), j2);
    line = [sPtY,sPtX,endPtY,endPtX];
    vecL = [line(3)-line(1),line(4)-line(2)];
    % get top and end pts
    if endPtY > extendBotPtVert(1)  % y is flipped in Matlab-- top is min
         extendBotPtVert(1) = endPtY;
         extendBotPtVert(2) = endPtX;
    end
    
    if endPtX > extendEndPtVert(2)
        extendEndPtVert(2) = endPtX;
        extendEndPtVert(1) = endPtY;
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

% scale back to the original image
extendTopPt = extendTopPt .* 2;
extendEndPt = extendEndPt .* 2;
extendBotPtVert = extendBotPtVert .*2;
extendEndPtVert = extendEndPtVert .*2;
figure; imshow(I); hold on;
eTop = [extendTopPt, extendBotPtVert];

plot([eTop(2),eTop(4)], [eTop(1),eTop(3)], 'r');
hold off;



