function preprocessROI(ROIfname)
close all;
I = imread(ROIfname);
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

% normalize values
Gmag = rescale(Gmag, 0, 1);
tmp = rescale(tmp, 0, 1);


% extract regions with Low Gmag and High intensity
Iupper = zeros(size(tmp));
Ilower = zeros(size(tmp));

gradThresUpperBone = 0.2;
valThresUpperBone = 0.7;

gradThresLowerBone = 0.3;
valThresLowerBone = 0.5;

uidx1 = Gmag(:,:) < gradThresUpperBone;
uidx2 = tmp(:,:) > valThresUpperBone; 

lidx1 = Gmag(:,:) > gradThresLowerBone;
lidx2 = tmp(:,:) > valThresLowerBone;

idx = uidx1 & uidx2;
% idx = uidx2;
Iupper(idx) = tmp(idx);
lidx = lidx1 & lidx2;
Ilower(lidx) = tmp(lidx);

% show the thresholded image and ask the user to draw a point
imshowpair(tmp, Gmag, 'montage');
bwUpper = imbinarize(Iupper);
bwUpper = bwmorph(bwUpper, 'bridge'); % bridge small gaps in the BW
bwUpper = bwmorph(bwUpper, 'spur');

CC = bwconncomp(bwUpper, 4) 
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
for i = 1:length(CC.PixelIdxList)
    pxList = CC.PixelIdxList{i};
    if isempty(find(pxList == ptIdx, 1))
        bwUpper(CC.PixelIdxList{i})= 0;
        continue;
    else
        disp('CC found!');
    end
end
% minLenThres = 50;
skel = bwskel(bwUpper);
X = bwmorph(skel, 'spur', inf);
imshow(X);
B = bwmorph(skel, 'branchpoints');
E = bwmorph(skel, 'endpoints');
[y,x] = find(E);
[y1,x1] = find(B);

figure; imshow(skel);hold on;
plot(x,y,'ro');
plot(x1,y1,'bo');
hold off;

B_loc = find(B);
E_loc = find(E);
% Dmask = false(size(skel));
D_list = cell(numel(x1),1);
% distToBranchPt = realmax;
for k = 1:numel(x1)
     D = bwdistgeodesic(skel,x1(k),y1(k));
     D_list{k} = D;
end

cutBranch = 0;
for k = 1:numel(x1)
    Dmask = false(size(skel));
    distToBranchPt = realmax;
%     D = bwdistgeodesic(skel,x1(k),y1(k));
    aE = 1; % find the associated end point of the branch
    for j = 1:numel(x)     
            isDirectConnect = true;
        for kk = 1:numel(x1)
            if kk == k
                continue;
            end
            
            if  D_list{k}(y(j),x(j)) > D_list{kk}(y(j),x(j))%make sure the branch pt does not go through another branch pt
                isDirectConnect = false;
                break;
            end
        end
        
        if isDirectConnect
            if  D_list{k}(y(j),x(j)) < distToBranchPt 
                distToBranchPt =  D_list{k}(y(j),x(j));
                aE = j;
            end
        end
    end
    aE
    % compute the property of the branch
    if y(aE)<y1(k) % if the branch is going up, remove it
        % make a horizontal cut
        cutBranch = aE;
        Dmask(min(y1(k),y(aE)):max(y1(k),y(aE)),:) = true;
        Bskel = skel & ~Dmask;
        Bskel = bwmorph(Bskel, 'spur');
        Bskel = bwmorph(Bskel, 'spur');
        Bskel = bwmorph(Bskel, 'spur');
        figure,imshow(Bskel);
        stats = regionprops('table',Bskel,'Centroid',...
            'MajorAxisLength','MinorAxisLength')
        break;
    end
    
%     Dmask(D_list{k}<distToBranchPt) = true;
%    distanceToBranchPt = min(D(E_loc));    
%    Dmask(D < distanceToBranchPt) =true;
end

% trace end points without cutBranch
skelInd = find(Bskel);
lastVecL = [0 1];
simpK = [];
% get end points
extendTopPt = [realmax,-realmax];
extendEndPt = [-realmax,-realmax];
%Bd = ind2sub(size(Bskel),skelInd);
for i = 1:numel(skelInd)-1
    j = skelInd(i); j2 = skelInd(i+1);
    [endPtY, endPtX] = ind2sub([size(Bskel,1),size(Bskel,2)], j);
    [sPtY,sPtX] = ind2sub([size(Bskel,1),size(Bskel,2)], j2);
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
[simpPtY,simpPtX] = ind2sub(size(Bskel), simpK);
plot(simpPtX,simpPtY, 'LineWidth', 3); 
plot([extendTopPt(2),extendEndPt(2)], [extendTopPt(1),extendEndPt(1)], 'r');
hold off;


% find the convex hull of points
% [k,av] = convhull(Bd);
% find the boundary (nonconvex) shape of points
% skelD = skel-Dmask;
% imshow(skelD);
hold all;

% bwUpper = bwUpper  bwUpperX;
imshow(bwUpper);


figure;
subplot(4, 1,1);
imshowpair(I, tmp, 'montage');
subplot(4,1,2);
imshowpair(tmp, Gmag, 'montage');
subplot(4,1,3);
% upper bone boundary
imshowpair(Gmag, Iupper, 'montage');
subplot(4,1,4);
% lower bone boundary
imshowpair(tmp, lidx1, 'montage');



% % find skeleton 
% SK = bwmorph(gpuArray(Iupper), 'skel', Inf);
% % SK = bwmorph(gpuArray(M), 'thin', Inf);
% 
% % find the contour of the image
% SK = bwperim(Iupper, 8);
% [B, L] = bwboundaries(Iupper);


