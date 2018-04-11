% function alpha = Fun_SmoothingMatting_CHU_V5(imageL, packL, imageR, packR)
% input: image trimap confidence alpha_hat
% 0722  try to getNewTrimap

% confidence=1-CalCost(input,F,B,output);
% conf_norm=(confidence-min(confidence(:)))/(max(confidence(:))-min(confidence(:)));

clear all
load('matlab data\test9_L2R_op')
load('matlab data\test9L_smooth_autotri.mat');
imageL = imread('test9L.png');
% imageL = IColor;  % original image
packL = pack;
load('matlab data\test9R_smooth_autotri.mat');
imageR = imread('test9R.png');
% imageR = IColor;
packR = pack; 
clear  pack  %IColor

% Parameters
omega    = 1e-1;
lambda   = 100;
[h,w,c] = size(imageR);
img_size = w * h;

  
% Normalize inputs
%  ****** left image *********
imageL      = double(imageL);  
alpha_hatL  = double(packL(:,:,1)); % baby2 size 415
confidenceL = double(packL(:,:,2));
alpha_hatL = alpha_hatL(1:h,1:w);     % baby2 different size
confidenceL = confidenceL(1:h,1:w);   % baby2 different size
trimapL = im2double(imread('test9L_autotri.png'));   % baby2 size 413
% trimapL = rgb2gray(trimapL);   % auto trimap does not need this line
clear packL
if max(imageL(:)) > 1
    imageL = imageL / 255;
end
if max(confidenceL(:)) > 1
    confidenceL = confidenceL / 255;
end
if max(alpha_hatL(:)) > 1
    alpha_hatL = alpha_hatL / 255;
end
consts_mapL = ~(trimapL(:,:,1) > 0 & trimapL(:,:,1) < max(trimapL(:))); 
            % 1 = foreground and background pixels, 0 = unknown pixels
%  ****** right image **********
imageR      = double(imageR);
alpha_hatR  = double(packR(:,:,1));
confidenceR = double(packR(:,:,2));
alpha_hatR = alpha_hatR(1:h,1:w);     % baby2 different size
confidenceR = confidenceR(1:h,1:w);   % baby2 different size
% trimapR     = packR(:,:,3);
% trimapR = trimapR(:,1:w2);
trimapR = im2double(imread('test9R_autotri.png')); 
% trimapR = rgb2gray(trimapR);   % auto trimap does not need this line
clear packR
if max(imageR(:)) > 1
    imageR = imageR / 255;
end
if max(confidenceR(:)) > 1
    confidenceR = confidenceR / 255;
end
if max(alpha_hatR(:)) > 1
    alpha_hatR = alpha_hatR / 255;
end
consts_mapR = ~(trimapR(:,:,1) > 0 & trimapR(:,:,1) < max(trimapR(:)));

%%  smoothing part
%  ***********  0725  update alpha_hat  *************
if vx(1,30)>0
   LAB = rgb2lab(imageR);
   Light = LAB(:,:,1);
   Light = Light/max(Light(:));
   a1 = alpha_hatR<1 & alpha_hatR>0;
   alp = alpha_hatR;
   alp(alp>0.4) = 1;   % threshold
   edg = edge(alp,'sobel');
   edg = imdilate(edg,ones(2));
   L1 = Light>0.45;     % threshold
   L1 = double(L1);
   L1(L1==1) = 5;
   a2 = a1-edg;
   a2(a2==1) = 5;
   a2(a2~=5) = 3;
   a3 = a2-L1;
   a3(a3~=0)=1;
   a3 =~a3;
   alpha_hatR(a3) = 1;
   L2 = Light<0.25;    % threshold
   L2(L2==1) = 5;
   a3 = a2-L2;
   a3(a3~=0)=1;
   a3 =~a3;
   alpha_hatR(a3) = max(alpha_hatR(a3)-0.2);  % threshold
else   
   LAB = rgb2lab(imageL);
   Light = LAB(:,:,1);
   Light = Light/max(Light(:));
   a1 = alpha_hatL<1 & alpha_hatL>0;
   alp = alpha_hatL;
   alp(alp>0.45) = 1;    % threshold
   edg = edge(alp,'sobel');
   edg = imdilate(edg,ones(2));  
   L1 = Light>0.5;     % threshold
   L1 = double(L1);
   L1(L1==1) = 5;
   a2 = a1-edg;
   a2(a2==1) = 5;
   a2(a2~=5) = 3;
   a3 = a2-L1;
   a3(a3~=0)=1;
   a3 =~a3;
   alpha_hatL(a3) = 1;
   L2 = Light<0.26;    % threshold
   L2(L2==1) = 5;
   a3 = a2-L2;
   a3(a3~=0)=1;
   a3 =~a3;
   alpha_hatL(a3) = max(alpha_hatL(a3)-0.3);  % threshold
end

% a1 = alpha_hatR<1 & alpha_hatR>0;
% alp = alpha_hatR;
% alp(alp>0.4) = 1;
% ed = edge(alp,'sobel');
% ed = imdilate(ed,ones(2));
% LAB = rgb2lab(imageR);
% Light = LAB(:,:,1);
% Light = Light/max(Light(:));
% L1 = Light>0.6;
% L1 = double(L1);
% L1(L1==1) = 5;
% a2 = a1-ed;
% a2(a2==1) = 5;
% a2(a2~=5) = 3;
% a3 = a2-L1;
% a3(a3~=0)=1;
% a3=~a3;
% alpha_hatR(a3) = 1;
% L2 = Light<0.16;
% L2(L2==1) = 5;
% a3 = a2-L2;
% a3(a3~=0)=1;
% a3=~a3;
% alpha_hatR(a3) = max(alpha_hatR(a3)-0.1);

% ******* decide left or right image first to get alpha values  *******
if vx(1,30)>0
    I1 = imageR;     % right image first, from right to left 
    I2 = imageL;  
%     [newtrimapR, alpha_hatR, confidenceR] = getNewTrimap_CHU(trimapR,alpha_hatR,alpha_hatL,confidenceR,vx);
%     consts = ~(newtrimapR(:,:,1) > 0 & newtrimapR(:,:,1) < max(newtrimapR(:)));
    trimap1 = newtrimapR;
%     [alpha_hatR, confidenceR] = getNewTrimap_CHU(trimapR,alpha_hatR,alpha_hatL,confidenceR,vx);
%     fR = confidenceR(:) .* ~consts_mapR(:);
     
    fR = confidenceR(:) .* ~consts(:);
    fL = confidenceL(:) .* ~consts_mapL(:);
    f = [fR; fL];
    alpha_hat = [alpha_hatR(:); alpha_hatL(:)];       
    consts_map = [consts(:);consts_mapL(:)];
%     consts_map = [consts_mapR(:);consts_mapL(:)];
%     consts = consts_mapR; 
else
    I1 = imageL;     % left image first  
    I2 = imageR; 
    [newtrimapL, alpha_hatL, confidenceL] = getNewTrimap_CHU(trimapL,alpha_hatL,alpha_hatR,confidenceL,vx);
    consts = ~(newtrimapL(:,:,1) > 0 & newtrimapL(:,:,1) < max(newtrimapL(:)));
    
%     [alpha_hatL, confidenceL] = getNewTrimap_CHU(trimapL,alpha_hatL,alpha_hatR,confidenceL,vx); 
%     fL = confidenceL(:) .* ~consts_mapL(:);
    trimap1 = newtrimapL;
    fL = confidenceL(:) .* ~consts(:);
    fR = confidenceR(:) .* ~consts_mapR(:);
    f = [fL; fR];    
    alpha_hat = [alpha_hatL(:); alpha_hatR(:)];    
    consts_map = [consts(:);consts_mapR(:)];        
%     consts_map = [consts_mapL(:);consts_mapR(:)]; 
%     consts = consts_mapL;
end

D = spdiags(consts_map, 0, 2*img_size, 2*img_size);
P = spdiags(f, 0, 2*img_size, 2*img_size);
K = lambda * D + omega * P;

% Generate matting Laplacian matrix
% L = getLaplacian_CHU_V5(I1, I2, consts,vx);
L = getLaplacian_CHU_V6(I1, I2,trimap1,edg,vx);

% Solve for alpha
x = (L + K) \ (K * alpha_hat(:));
alpha = max(min(reshape(x,h,2*w),1),0);
alpha1 = alpha(:,1:w);
k = alpha1<0.1 ;
alpha1(k)=0;
k = alpha1>0.9;
alpha1(k) = 1;
figure,imshow(alpha)
% figure,imshow(alpha1)

WhiteBack = ones(h,w,c);
ComF = bsxfun(@times,I1,alpha1)+bsxfun(@times,(1-alpha1),WhiteBack);
% ComF = ComF/255;
% figure,imshow(ComF)

