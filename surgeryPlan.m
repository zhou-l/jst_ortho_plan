% stretchAlpha is the angle (in degrees) of stretching
function Iplan = surgeryPlan( xRayInput, xRayUpperFname, xRayLowerFname, segMaskUpperBWFile, segMaskLowerBWFile, stretchAlpha, isRightInImg)
if nargin < 6
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
Ifull = imread(xRayInput);
hI = imshow(Ifull);hold on;
% alpha(hI, 0.7);
hold on;
isLower = false;
[hAxUpper, vAxUpper, boundUpper, kUpper] = findBoneAxis(xRayUpperFname, segMaskUpperBWFile, isLower);
hold on;
isLower = true;
[hAxLower, vAxLower, boundLower, kLower] = findBoneAxis(xRayLowerFname, segMaskLowerBWFile, isLower);

scale = 600;

% plot required axes
% upper 
plot([hAxUpper(:,1)-scale .* hAxUpper(:,3), hAxUpper(:,1)+scale .* hAxUpper(:,3)],... 
     [hAxUpper(:,2)-scale .* hAxUpper(:,4), hAxUpper(:,2)+scale .* hAxUpper(:,4)], 'LineWidth', 2);
plot([vAxUpper(:,1)-scale .* vAxUpper(:,3), vAxUpper(:,1)+scale .* vAxUpper(:,3)],...
     [vAxUpper(:,2)-scale .* vAxUpper(:,4), vAxUpper(:,2)+scale .* vAxUpper(:,4)], 'LineWidth', 2);
% lower
plot([hAxLower(:,1)-scale .* hAxLower(:,3), hAxLower(:,1)+scale .* hAxLower(:,3)],...
     [hAxLower(:,2)-scale .* hAxLower(:,4), hAxLower(:,2)+scale .* hAxLower(:,4)], 'LineWidth', 2);
plot([vAxLower(:,1)-scale .* vAxLower(:,3), vAxLower(:,1)+scale .* vAxLower(:,3)],...
     [vAxLower(:,2)-scale .* vAxLower(:,4), vAxLower(:,2)+scale .* vAxLower(:,4)], 'LineWidth', 2);
 
% find the intersection of two vertical axes
dx = vAxLower(1) - vAxUpper(1);
dy = vAxLower(2) - vAxUpper(2);
det = vAxLower(3) * vAxUpper(4) - vAxLower(4) * vAxUpper(3);
if abs(det) < 1e-8
    warning('This is a rare case that two axes do not intersect! Quit!');
    return ;
end

% find the intersection point
u = (dy * vAxLower(3) - dx * vAxLower(4)) / det; 
v = (dy * vAxUpper(3) - dx * vAxUpper(4)) / det;

intPt = [vAxLower(1) + v * vAxLower(3),vAxLower(2) + v * vAxLower(4)];
% draw horizontal axis at this point
hPt0 = [intPt(1)-scale* hAxUpper(:,3), intPt(2)-scale* hAxUpper(:,4)];
hPt1 = [intPt(1)+scale* hAxUpper(:,3), intPt(2)+scale* hAxUpper(:,4)];
plot([hPt0(1),hPt1(1)], [hPt0(2), hPt1(2)], 'LineWidth', 3);

% find intersection of the horizontal line and the boundary
baselineIntPts = [];
for i = 1:length(kUpper)-1
    bPt0 = [boundUpper(kUpper(i),2), boundUpper(kUpper(i),1)];
    bPt1 = [boundUpper(kUpper(i+1),2), boundUpper(kUpper(i+1),1)];
    if intersectLL2D(hPt0, hPt1, bPt0, bPt1)
        lines = [hPt0; hPt1; bPt0; bPt1];
        pt = linlinintersect(lines);
        baselineIntPts(end+1,:) = pt;
        plot(pt(1), pt(2), 'o', 'LineWidth', 2);
    end
end

if isempty(baselineIntPts)
    disp('Cannot find cut intersection point!');
    return;
end
% draw cut path
cutEntryPt = baselineIntPts(1,:);
if size(baselineIntPts,1) == 2
    if baselineIntPts(1,1) > baselineIntPts(2,1)  % element 1 is on the right 
        if zIsRightInImg
            cutEntryPt = baselineIntPts(1,:);
            cutExitPt = baselineIntPts(2,:);
        else
            cutEntryPt = baselineIntPts(2,:);
            cutExitPt = baselineIntPts(1,:);
        end
    else
        if zIsRightInImg
            cutEntryPt = baselineIntPts(2,:);
            cutExitPt = baselineIntPts(1,:);
        else
            cutEntryPt = baselineIntPts(1,:);
            cutExitPt = baselineIntPts(2,:);
        end
    end
end
% The cutting ray: 
if zIsRightInImg
%     Rcy = cutEntryPt + t * sin(zalpha/180*pi);
%     Rcx = cutEntryPt - t * cos(zalpha/180*pi);
      Rc = [cutEntryPt(1) cutEntryPt(2) -cos(zStretchAlpha/180*pi) -sin(zStretchAlpha/180*pi) ];
else
%     Rcy = cutEntryPt + t * sin(zalpha/180*pi);
%     Rcx = cutEntryPt + t * cos(zalpha/180*pi);  
    Rc = [cutEntryPt(1) cutEntryPt(2)  cos(zStretchAlpha/180*pi) -sin(zStretchAlpha/180*pi)];
end
% draw the cutting ray
cPt0 = [cutEntryPt(1)-scale* Rc(:,3), cutEntryPt(2)-scale* Rc(:,4)];
cPt1 = [cutEntryPt(1)+scale* Rc(:,3), cutEntryPt(2)+scale* Rc(:,4)];
plot([cPt0(1),cPt1(1)], [cPt0(2), cPt1(2)], 'LineWidth', 3);

cutTurnPt = cutExitPt;
for i = 1:length(kUpper)-1
    bPt0 = [boundUpper(kUpper(i),2), boundUpper(kUpper(i),1)];
    bPt1 = [boundUpper(kUpper(i+1),2), boundUpper(kUpper(i+1),1)];
    if intersectLL2D(cPt0, cPt1, bPt0, bPt1)
        lines = [cPt0; cPt1; bPt0; bPt1];
        pt = linlinintersect(lines);
        
        if norm(pt - cutEntryPt) < 1e-4
            continue;
        else
            cutTurnPt = pt;
            plot(pt(1), pt(2), 'o', 'LineWidth', 2);
        end
    end
end

% connect the cut path: entry -> turn -> exit pts
plot([cutEntryPt(1) cutTurnPt(1) cutExitPt(1)], [cutEntryPt(2) cutTurnPt(2) cutExitPt(2)], 'r', 'LineWidth', 3);


