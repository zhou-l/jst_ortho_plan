function eBottom = preprocessROIbottom(ROIfname, zisRight)
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
% corners = detectHarrisFeatures(tmp);
% imshow(tmp); hold on;
% plot(corners.selectStrongest(50));

BWE = edge(tmp, 'Canny');
imshowpair(tmp, BWE, 'montage');

% for the bottom part, Canny edge detection does a fairly good job! 
% We could start from there directly
bwLower = bwmorph(BWE, 'clean');
bwLower = bwmorph(bwLower, 'spur');
 figure;imshow(bwLower); 
CC = bwconncomp(bwLower) 

skel = bwLower;
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

% and search for the vertical CC
vertCC = 0;
minLenThres = 15;
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
        skel(CC.PixelIdxList{i}) = true;
        vertCC = i;
    end 
end
% skel = bwmorph(skel, 'thicken');
% skel = bwmorph(skel,'bridge');
% 
%  E = bwmorph(skel, 'endpoints');
%  imshow(skel); hold on;
% [y1,x1] = find(E);
% plot(x1,y1,'ro');

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

% scale back to the original image

extendTopPt = extendTopPt .*2;
extendTopPtVert = extendTopPtVert .*2;
figure; imshow(I); hold on;
eBottom = [extendTopPt extendTopPtVert];
plot([extendTopPt(2),extendTopPtVert(2)], [extendTopPt(1),extendTopPtVert(1)], 'r');
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


