% function alpha = Fun_SmoothingMatting_CHU_V5(imageL, packL, imageR, packR)
% input: image trimap confidence alpha_hat
% 0727 update alpha before smooth 

% confidence=1-CalCost(input,F,B,output);
% conf_norm=(confidence-min(confidence(:)))/(max(confidence(:))-min(confidence(:)));

clear all
addpath('0925')
load('matlab data\test20_R2L_op')
load('matlab data\test20L_lb.mat');
imageL = im2double(imread('test20L.png'));
trimapL = im2double(imread('test20L_tri.png')); 
% packL = pack;
alpha_hatL = alpha;
confidenceL = confidence;
%alpha_hatL = im2double(imread('test20L_cs.png'));
load('matlab data\test20R_lb.mat');
imageR = im2double(imread('test20R.png'));
trimapR = im2double(imread('test20R_tri.png'));
%alpha_hatR = im2double(imread('test20R_cs.png'));
% packR = pack; 
alpha_hatR = alpha;
confidenceR = confidence;
% clear  pack  
% Parameters
omega    = 1e-1;
lambda   = 100;
[h,w,c] = size(imageR);
img_size = w * h;

  
% Normalize inputs
%  ****** left image *********
% imageL      = double(imageL);  
% alpha_hatL  = double(packL(:,:,1)); % baby2 size 415
% confidenceL = double(packL(:,:,2));
alpha_hatL = alpha_hatL(1:h,1:w);     % baby2 different size
confidenceL = confidenceL(1:h,1:w);   % baby2 different size
% trimapL = im2double(imread('test9L_autotri.png'));   % baby2 size 413
% trimapL = rgb2gray(trimapL);   % auto trimap does not need this line
% clear packL
if max(imageL(:)) > 1
    imageL = imageL / 255;
end
if max(confidenceL(:)) > 1
    confidenceL = confidenceL / 255;
end
% if max(alpha_hatL(:)) > 1
%     alpha_hatL = alpha_hatL / 255;
% end
consts_mapL = ~(trimapL(:,:,1) > 0 & trimapL(:,:,1) < max(trimapL(:))); 
            % 1 = foreground and background pixels, 0 = unknown pixels
%  ****** right image **********
% imageR      = double(imageR);
% alpha_hatR  = double(packR(:,:,1));
% confidenceR = double(packR(:,:,2));
alpha_hatR = alpha_hatR(1:h,1:w);     % baby2 different size
confidenceR = confidenceR(1:h,1:w);   % baby2 different size
% trimapR     = packR(:,:,3);
% trimapR = trimapR(:,1:w2);
% trimapR = im2double(imread('test9R_autotri.png')); 
% trimapR = rgb2gray(trimapR);   % auto trimap does not need this line
% clear packR
if max(imageR(:)) > 1
    imageR = imageR / 255;
end
if max(confidenceR(:)) > 1
    confidenceR = confidenceR / 255;
end
% if max(alpha_hatR(:)) > 1
%     alpha_hatR = alpha_hatR / 255;
% end
trimapR(alpha_hatR>0.96)=1;
consts_mapR = ~(trimapR(:,:,1) > 0 & trimapR(:,:,1) < max(trimapR(:)));

%%  smoothing part
%  ***********  0725  update alpha_hat  *************
if vx(1,30)>0
     figure,imshow(alpha_hatR)
   LAB = rgb2lab(imageR);
   Light = LAB(:,:,1);
   Light = Light/max(Light(:));
   a1 = alpha_hatR<0.98 & alpha_hatR>0.5;
   alp = alpha_hatR;
% %    alp(alp>0.6) = 1;   % threshold
   edg = edge(alp,'sobel');
    edg = imdilate(edg,ones(3));
   L1 = Light>0.4; %|( Light>0.18 & Light<0.5) ;     % threshold
   L1 = double(L1);
   L1(L1==1) = 5;
   a2 = a1-edg;
   a2(a2==1) = 5;
   a2(a2~=5) = 3;
   a3 = a2-L1;
   a3(a3~=0)=1;
   a3 =~a3;
   alpha_hatR(a3) = 1;
  
% %    L2 = Light<0.1;    % threshold
% %    L2(L2==1) = 5;
% %    a3 = a2-L2;
% %    a3(a3~=0)=1;
% %    a3 =~a3;
% %    alpha_hatR(a3) = max(alpha_hatR(a3)-0.1,0);  % threshold
else  
    figure, imshow(alpha_hatL)
   LAB = rgb2lab(imageL);
   Light = LAB(:,:,1);
   Light = Light/max(Light(:));
   a1 = alpha_hatL<0.98 & alpha_hatL>0.5;
   alp = alpha_hatL;
% %    alp(alp>0.6) = 1;    % threshold
   edg = edge(alp,'sobel');
   edg = imdilate(edg,ones(3));  
   L1 = Light>0.2;     % threshold
   L1 = double(L1);
   L1(L1==1) = 5;
   a2 = a1-edg;
   a2(a2==1) = 5;
   a2(a2~=5) = 3;
   a3 = a2-L1;
   a3(a3~=0)=1;
   a3 =~a3;
   alpha_hatL(a3) = 1;
% %    L2 = Light<0.1;    % threshold
% %    L2(L2==1) = 5;
% %    a3 = a2-L2;
% %    a3(a3~=0)=1;
% %    a3 =~a3;
% %    alpha_hatL(a3) = max(alpha_hatL(a3)-0.1,0);  % threshold
end

% ******* decide left or right image first to get alpha values  *******
fL = confidenceL(:) .* ~consts_mapL(:);
fR = confidenceR(:) .* ~consts_mapR(:);
if vx(1,30)>0
    I1 = imageR;     % right image first, from right to left 
    I2 = imageL;       
    f = [fR; fL];
    alpha_hat = [alpha_hatR(:); alpha_hatL(:)];       
    consts_map = [consts_mapR(:);consts_mapL(:)];
    consts = consts_mapR;     % when use original trimap, give up imerode 
    trimap1 = trimapR;
else
    I1 = imageL;     % left image first  
    I2 = imageR;     
    f = [fL; fR];    
    alpha_hat = [alpha_hatL(:); alpha_hatR(:)];    
    consts_map = [consts_mapL(:);consts_mapR(:)];  
    consts = consts_mapL;   
    trimap1 = trimapL;
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
% k = alpha1<0.1 ;
% alpha1(k)=0;
% k = alpha1>0.9;
% alpha1(k) = 1;
% figure,imshow(alpha)
WhiteBack = ones(h,w,c);
ComF = bsxfun(@times,I1,alpha1)+bsxfun(@times,(1-alpha1),WhiteBack);
% ComF = ComF/255;
 figure,imshow(ComF)
 figure,imshow(alpha)
% imwrite(alpha1,'test2L_knn_up.png','png')