function [flow,warpmap] = opticalflow(im1,im2,par)
%addpath('mex');
addpath('OpticalFlow\mex');
im1 = im2double(im1);
im2 = im2double(im2);


% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
% par.alpha = 0.012;
% par.ratio = 0.75;
% par.minWidth = 20;
% par.nOuterFPIterations = 7;
% par.nInnerFPIterations = 1;
% par.nSORIterations = 30;

para = [par.alpha,par.ratio,par.minWidth,par.nOuterFPIterations,par.nInnerFPIterations,par.nSORIterations];

% this is the core part of calling the mexed dll file for computing optical flow
% it also returns the time that is needed for two-frame estimation
tic;
[vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);%warpI2 is warp im2 to im1.Im1 is the main image.
toc


%figure;imshow(im1);
figure;imshow(warpI2);title('warp image');

% visualize flow field
flow(:,:,1) = vx;
flow(:,:,2) = vy;
imflow = flowToColor(flow);

% flow1(:,:,1) = vx1;
% flow1(:,:,2) = vy1;
% imflow1 = flowToColor(flow1);

figure;imshow(imflow);title('flow image')
warpmap=im2uint8(warpI2);
end
% figure;imshow(imflow1);
