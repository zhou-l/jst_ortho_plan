% Find the axes of the bone:
% The output are stored as lines: [x0 y0 v0.x v0.y]
% horAxis  - horizontal axis
% vertAxis - vertical axis
function [horAxis, vertAxis, Bd, simpK] = findBoneAxis(xRayFname, segMaskBWFile, isLower)

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
% SK = bwmorph(gpuArray(M), 'thin', Inf);

% find the contour of the image
SK = bwperim(M, 8);
[B, L] = bwboundaries(M);
% SK = bwskel(M);
% overlay with label 
% imshow(labeloverlay(I, SK, 'Transparency', 0));

% overlay with transparency

% find the largest object
objI = -1;
maxBoundPxCnt = -1;
for k = 1:length(B)
    if maxBoundPxCnt < length(B{k})
        maxBoundPxCnt = length(B{k});
        objI = k;
    end
end

hI = imshow(I);hold on;
alpha(hI, 0.3);
Bd = B{objI};
% plot(Bd(:,2), Bd(:,1), 'w', 'LineWidth', 1);

% find the convex hull of points
% [k,av] = convhull(Bd);
% find the boundary (nonconvex) shape of points
[k,av] = boundary(Bd, 0.9);
% plot(Bd(k,2),Bd(k,1), 'LineWidth', 4);
% pgon = polyshape(Bd(k,2), Bd(k,1));

% find edges that are horizontal...
maxHorLlen = -1;

% do smoothing/simplification
lastVecL = [0 1];
simpK = [];
for i = 1:length(k)-1
    j = k(i); j2 = k(i+1);
    line = [Bd(j,2),Bd(j,1),Bd(j2,2),Bd(j2,1)];
    vecL = [line(3)-line(1),line(4)-line(2)];
    if norm(vecL) > 1e-4 
        vecL = vecL / norm(vecL);
    else
        vecL = [0 0];
    end
    % angle difference is small
    if dot(lastVecL,vecL)>cos(pi/30)
        continue;     
    else
        simpK(end+1,:) = j;
        lastVecL = vecL;     
    end
end

plot(Bd(simpK,2),Bd(simpK,1), 'LineWidth', 3);
for i = 1:length(simpK)-1
    j = simpK(i); j2 = simpK(i+1);
    line = [Bd(j,2),Bd(j,1),Bd(j2,2),Bd(j2,1)];
    vecL = [line(3)-line(1),line(4)-line(2)];
    normL = norm(vecL);
    normVecL = vecL ./ normL;
    % almost parallel to the horizontal vector
    if isLower 
        if dot(vecL/norm(vecL), [-1,0])< cos(pi/20)
            continue;
        end
    else
        if dot(vecL, [1,0])/norm(vecL) < cos(pi/20)
            continue;
        end
    end
 
    midPt = [(line(3)+line(1))/2, (line(4)+line(2))/2];
    % the normal vector!
    pVL = [vecL(2),-vecL(1)];
 
%     pVL = pVL ./ norm(pVL);

    x = int16([midPt(1),midPt(1)+pVL(1)*11]);
    y = int16([midPt(2),midPt(2)+pVL(2)*11]);
    
      plot(x, y);
    if maxHorLlen < normL
        maxHorLlen = normL;
        horAxis = [midPt(1) midPt(2) normVecL(1) normVecL(2)];
        normPVL = norm(pVL);
        vertAxis = [midPt(1) midPt(2) pVL(1)./normPVL pVL(2)./normPVL];
    end
end

plot([vertAxis(:,1), vertAxis(:,1)+100 .* vertAxis(:,3)], [vertAxis(:,2), vertAxis(:,2)+100 .* vertAxis(:,4)], 'LineWidth', 2);
plot([horAxis(:,1), horAxis(:,1)+100 .* horAxis(:,3)], [horAxis(:,2), horAxis(:,2)+100 .* horAxis(:,4)], 'LineWidth', 2);
