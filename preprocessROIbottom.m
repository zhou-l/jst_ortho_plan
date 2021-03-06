function eBottom = preprocessROIbottom(ROIfname, zisRight)
global npix10cm;
if isempty(npix10cm)
    znipx10cm = 724; % a default value
else
    znipx10cm = npix10cm;
end
%% find the horizontal line of the bottom bone at the joint
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

% tmp = imdiffusefilt(tmp);
[Gmag, Gdir]=imgradient(tmp);

Gmag = imdiffusefilt(Gmag);
imshowpair(tmp, Gmag, 'montage');

Gmag = rescale(Gmag, 0, 1);
tmp = rescale(tmp, 0, 1);
% % try detecting corner points...
%% extract high intensity boundary of the bottom  part
gradThresBone = 0.2;
valThresBone = 0.7;

uidx1 = Gmag(:,:) > gradThresBone;
uidx2 = tmp(:,:) > valThresBone; 
idx = uidx2; %uidx1 & uidx2; 
tmpX = zeros(size(tmp));
tmpX(idx) = tmp(idx);

[BWE,threshout] = edge(tmp, 'Canny');
[BWE,threshout2] = edge(tmp, 'Canny', [threshout(2)*0.6,threshout(2)*1.5]);
figure; imshowpair(tmpX, BWE, 'montage');
BWEo = BWE & imbinarize(tmpX);
figure; imshowpair(BWEo, BWE, 'montage');

% for the bottom part, Canny edge detection does a fairly good job! 
% We could start from there directly
bwLower = bwmorph(BWE, 'clean');
bwLower = bwmorph(bwLower, 'spur');

% bwLower = bwmorph(bwLower, 'bridge');
 figure;imshow(bwLower); 
CC = bwconncomp(bwLower) 

skel = false(size(tmp));
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
mainCC = idx;

% % Ask the user to choose the bottom boundary (mainCC)
% % h = drawpoint;
%   [xi, yi, but] = ginput(1);
%   xi = int16(xi);
%   yi = int16(yi);
%   ptIdx = int16((xi-1)*size(tmp,1)+yi);
%   ptIdx2 = sub2ind(size(tmp), yi, xi);
%   
% mainCC = 0;
% for i = 1:length(CC.PixelIdxList)
%     pxList = CC.PixelIdxList{i};
%     if isempty(find(pxList == ptIdx, 1))
%         skel(CC.PixelIdxList{i})= 0;
%         continue;
%     else
%         disp('CC found!');
%         mainCC = i;
%     end
% end
if mainCC < 1
    warning('Failed to find the bottom boundary!');
    return;
end
skel(CC.PixelIdxList{mainCC}) = true;

% and search for the vertical CC
vertCC = 0;
minLenThres = 15;
strongestEGmag = -realmax;
for i = 1:length(CC.PixelIdxList)
    if mainCC == i
        continue;
    end
    pxList = CC.PixelIdxList{i};
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
    if numel(pxList) <= minLenThres
        continue;
    end
    
    if isRight
        theta = acos(dot(dir,[1,0]));
    else
        theta = acos(dot(dir,[-1,0]));
    end
    if theta - pi/2 < 0 && -theta + pi/2 < pi/4 % The vertical edge should go from bottom right to top left
        eGmag = mean(Gmag(CC.PixelIdxList{i}));
        if eGmag > strongestEGmag % and the edge should be strong
           strongestEGmag = eGmag;
           vertCC = i;
       end
    end 
end
skel(CC.PixelIdxList{vertCC}) = true;
figure; imshow(skel);

% skel = bwmorph(skel, 'thicken');
% skel = bwmorph(skel,'bridge');
% 
%  E = bwmorph(skel, 'endpoints');
%  imshow(skel); hold on;
% [y1,x1] = find(E);
% plot(x1,y1,'ro');

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


lastVecL = [0 1];
simpK = [];
% get end points
extendTopPtVert = [realmax,-realmax];
extendEndPtVert = [-realmax,-realmax];
%Bd = ind2sub(size(Bskel),skelInd);
skelInd = CC.PixelIdxList{vertCC};
for i = 1:numel(skelInd)-1
    j = skelInd(i); j2 = skelInd(i+1);
    [endPtY, endPtX] = ind2sub(size(tmp), j);
    [sPtY,sPtX] = ind2sub(size(tmp), j2);
    line = [sPtY,sPtX,endPtY,endPtX];
    vecL = [line(3)-line(1),line(4)-line(2)];
    % get top and end pts
    if endPtY < extendTopPtVert(1)  % y is flipped in Matlab-- top is min
         extendTopPtVert(1) = endPtY;
         extendTopPtVert(2) = endPtX;
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

%% drop the line a bit to find the width of the bottom bone
eBottom2 = [extendTopPt extendTopPtVert];
dropDist = znipx10cm * 0.05 / 2; % drop down by 5mm--- and consider the scaling factor of 1/2
% find the perpendicular line to eBottom
dBot = normalize(eBottom2(3:4)-eBottom2(1:2),'norm',2);
if ~isRight
    dBot = -dBot;
end
dVert = [dBot(2), -dBot(1)];
dropPt = dVert .* dropDist + 0.5 .* (eBottom2(1:2) + eBottom2(3:4));
dropLs = dropPt - dBot * 50;
dropLt = dropPt + dBot * 50;
% find the intersections of the line and the bone 
% do ray marching for the outer side
stepSize = 1;
s = dropLt;
eBottom = [dropLs dropLt];
outerVertBoundary = CC.PixelIdxList{vertCC};
nstep = 500;
for i = 1:nstep
    if s(1) > size(tmp,1) || s(1) < 1 || s(2) > size(tmp,2) || s(2) < 1
        break;
    end
    s = s + dBot * stepSize;
    idxs = sub2ind(size(tmp), int16(s(1)), int16(s(2)));
    hitPts = find(outerVertBoundary==idxs, 1);
    if ~isempty(hitPts)
        eBottom(3:4) = int16(s);
        break;
    end
end
% ray marching for the inner side
s = dropLt;
innerBoundary = CC.PixelIdxList{mainCC};
for i = 1:nstep
    if s(1) > size(tmp,1) || s(1) < 1 || s(2) > size(tmp,2) || s(2) < 1
        break;
    end
    s = s - dBot * stepSize;
    idxs = sub2ind(size(tmp), int16(s(1)), int16(s(2)));
    hitPts = find(innerBoundary==idxs, 1);
    if ~isempty(hitPts)
        eBottom(1:2) = int16(s);
        break;
    end
end
figure; imshow(tmp); hold on;
plot([extendTopPt(2),extendTopPtVert(2)], [extendTopPt(1),extendTopPtVert(1)], 'b');
plot([eBottom(2),eBottom(4)], [eBottom(1),eBottom(3)],'r');
hold off;

% scale back to the original image

figure; imshow(I); hold on;
eBottom = 2 .* eBottom;
eBottom2 = 2 .* eBottom2;
plot([eBottom2(2),eBottom2(4)], [eBottom2(1),eBottom2(3)], 'b');
plot([eBottom(2),eBottom(4)], [eBottom(1),eBottom(3)],'r');
hold off;

% figure;
% subplot(4, 1,1);
% imshowpair(I, tmp, 'montage');
% subplot(4,1,2);
% imshowpair(tmp, Gmag, 'montage');
% subplot(4,1,3);
% % upper bone boundary
% imshowpair(Gmag, Iupper, 'montage');
% subplot(4,1,4);
% % lower bone boundary
% imshowpair(tmp, lidx1, 'montage');


