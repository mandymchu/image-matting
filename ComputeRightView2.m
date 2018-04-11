%  compute the right eye view from left view and disparity map

clear all
close all

I = im2double(imread('baby2_R_tri_ori.png')); 
% Disp = double(imread('baby3_R.png')); 
[IH , IW,ICmode]= size(I) ;
% load baby2_LtoR_vx;   % load optical flow result 
% vxL = vx;   vyL = vy;
load baby2_RtoL_vx;      % load baby2 optical flow result
% vxR = vx;   vyR = vy;
% vxerr = vxL+vxR;   vyerr = vyL+vyR;     %  error  
% % vxL = vxL-vxerr/2;  vyL = vyL-vyerr/2;  
% vx = vxR-vxerr/2;  vy = vyR-vyerr/2;

R1 = I(:,:,1);
G1 = I(:,:,2);
B1 = I(:,:,3);
% R2 = zeros(IH,IW);
R2 = ones(IH,IW);
R2 = (192/255)*R2;
G2 = R2;
B2 = R2;

% for i = 1:IH
%    for j = 1:IW
%        d = Disp(i,j)/4;
%        d = round(d);
%        t = j-d;
%        if t>0
%          R2(i,t) = R1(i,j);
%          G2(i,t) = G1(i,j);
%          B2(i,t) = B1(i,j);
%        end
%    end
% end


for i = 1:IH
    for j = 1:IW
     dx = round(vx(i,j));
     dy = vy(i,j);
     dy = round(dy);
      m = i+dy;
      n = j+dx;
     if m>0&&n>0&&m<(IH+1)&&n<(IW+1)
        R2(m,n) = R1(i,j); 
        G2(m,n) = G1(i,j);
        B2(m,n) = B1(i,j); 
     end
     
%        if y>0
%          R2(i,y) = R1(i,j);
%          G2(i,y) = G1(i,j);
%          B2(i,y) = B1(i,j);
%        end      
    end 
end

I2(:,:,1) = R2;
I2(:,:,2) = G2;
I2(:,:,3) = B2;
figure,imshow(I2)  
% imwrite(I2,'feather36_tri.png','png');

% load the original right eye view
I3 = im2double(imread('baby2_L.png'));
for i = 1:3
   I4(:,:,i) = I2(:,:,i).*I3(:,:,i); 
end
% subplot(3,1,1)
% imshow(I3)
% subplot(3,1,2)
% imshow(I2)
% subplot(3,1,3)
% imshow(I4)
figure,imshow(I4)

%%  

% I1 = im2double(imread('baby2_L.png'));   
% I2 = im2double(imread('baby2_R.png')); 
% I3 = double(imread('disp1.png'));   % left disparity map
% I4 = double(imread('disp5.png'));   % right disparity map
% I5 = im2double(imread('baby2_L_tri.png'));  % left trimap
% [IH , IW,~]= size(I1) ;
% 
% % R1 = I5(:,:,1);
% % G1 = I5(:,:,2);
% % B1 = I5(:,:,3);
% % R2 = ones(IH,IW);
% % R2 = (192/255)*R2;
% % G2 = R2;
% % B2 = R2;
% 
% I8 = im2double(imread('baby2_L_Alpha.png'));
% I9 = im2double(imread('baby2_R_Alpha_ori.png'));
% [IH2, IW2] = size(I8);
% I8new = I8;
% I9new = I9;
% 
% for i = 2:2:IH
%    for j = 2:2:IW-1
%         i2 = i/2;
%         j2 = j/2;
%        d = I3(i,j)/2.7;
%        d2 = d/2;
%        d = round(d);
%        d2 = round(d2);
%        t = j-d; 
%        t2 = j2-d2;
%        if t>0&&t2>0
%            m = I3(i,j);
%            n = I4(i,t);  
%            if m==n
%               I8new(i2,j2) = max(I8(i2,j2),I9(i2,t2));
%               I9new(i2,t2) = max(I8(i2,j2),I9(i2,t2));
%            else
%            I8new(i2,j2) = I8(i2,j2)*(m>n)+I9(i2,t2)*(n>m);
%            I9new(i2,t2) = I8(i2,j2)*(m>n)+I9(i2,t2)*(n>m);
%            end
% %          R2(i,t) = R1(i,j);
% %          G2(i,t) = G1(i,j);
% %          B2(i,t) = B1(i,j);
%         
%        end
%    end
% end
% figure,imshow(I8);
% figure,imshow(I9);
% figure,imshow(I8new);
% figure,imshow(I9new);

% I6(:,:,1) = R2;  % right trimap
% I6(:,:,2) = G2;
% I6(:,:,3) = B2;
% figure,imshow(I6); 
% 
% for i = 1:3
%    I7(:,:,i) = I2(:,:,i).*I6(:,:,i); 
% end
% figure,imshow(I7);



