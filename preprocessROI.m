function [eTop,eBot] = preprocessROI(ROIfname, zisRight)
close all;
if nargin < 2
    isRight = true;
else
    isRight = zisRight;
end
eTop = preprocessROItop2(ROIfname, isRight);
eBot = preprocessROIbottom(ROIfname, isRight);

I = imread(ROIfname);
figure;imshow(I);
hold on;
plot([eTop(2),eTop(4)],[eTop(1),eTop(3)],'r');
plot([eBot(2),eBot(4)],[eBot(1),eBot(3)],'b');
hold off;


