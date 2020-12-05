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
Dmask = false(size(skel));
for k = 1:numel(x1)
    D = bwdistgeodesic(skel,x1(k),y1(k));
    distanceToBranchPt = min(D(B_loc));
    Dmask(D < distanceToBranchPt) =true;
end
skelD = skel-Dmask;
imshow(skelD);
hold all;


imshow(out1)
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


