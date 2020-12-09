% stretchAlpha is the angle (in degrees) of stretching
function Iplan = surgeryPlanROI( xRayInputfName )
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
%% 1. ask the user to draw a roi rectangle for the 10cm measurement
hI = imshow(Ifull);
roi10cm = drawrectangle;
pos10cm = customWait(roi10cm); % pos is [x, y, width, height]
npix10cm = pos10cm(4); % get the length in pixels of 10cm

%% 2. ask the user to interactively draw a roi rectangle to select the joint of interest

roi = drawrectangle;
pos = customWait(roi);

% % crop the ROI image
% Ir = zeros(pos(4),pos(3));
Ir = Ifull(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3));
imwrite(Ir,'userROI.png');
[eTop, eBot] = preprocessROI('userROI.png');% line segments are stored as [ys xs yt xt]

% transform the lines back to the full image
eTop(1) = eTop(1) + pos(2); eTop(3) = eTop(3) + pos(2); 
eTop(2) = eTop(2) + pos(1); eTop(4) = eTop(4) + pos(1); 

eBot(1) = eBot(1) + pos(2); eBot(3) = eBot(3) + pos(2); 
eBot(2) = eBot(2) + pos(1); eBot(4) = eBot(4) + pos(1); 

imshow(Ifull);hold on;
plot([eTop(2),eTop(4)],[eTop(1),eTop(3)],'r');
plot([eBot(2),eBot(4)],[eBot(1),eBot(3)],'b');

%% 3. find the top vertical line
%% 3.1 find the 8cm on top of the mid-point of eTop
eTopCtr = [0.5 * (eTop(1) + eTop(3)), 0.5 * (eTop(2) + eTop(4))];
eTop8cm = int16(-0.8 * npix10cm + eTopCtr(1));
eTop13cm = int16(-1.3 * npix10cm + eTopCtr(1));

Itop = Ifull(eTop13cm:eTop8cm,int16(0.5*size(Ifull,2)):size(Ifull,2));
figure, imshow(Itop);
imwrite(Itop,'topBoneVertLineROI.png');
eTopVert = processTopBoneVertLine('topBoneVertLineROI.png');

% eBot(2) = eBot(2) + pos(1); eBot(4) = eBot(4) + pos(1); 

imshow(Ifull);hold on;
plot([eTop(2),eTop(4)],[eTop(1),eTop(3)],'r');
plot([eBot(2),eBot(4)],[eBot(1),eBot(3)],'b');
plot([eBot(2),eBot(4)],[eBot(1),eBot(3)],'b');
%%
%%4. find the bottom vertical line

I_xRayInput = Ifull;
[Gmag, Gdir] = imgradient(Ifull,'sobel');
imshow(Gmag);
% alpha(hI, 0.7);
hold on;
isLower = false;
[hAxUpper, vAxUpper, boundUpper, kUpper] = findUpperBoneAxis(I_xRayInput, I_upperBoneROI);
hold on;
isLower = true;
[hAxLower, vAxLower, boundLower, kLower] = findLowerBoneAxis(I_xRayInput, I_lowerBoneROI);
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