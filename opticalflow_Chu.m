clear all
close all
IL = im2double(imread('baby2_L_Alpha.png'));  % left view
IR = im2double(imread('baby2_R_Alpha_ori.png'));  % right view
[IH , IW]= size(IL); 
load baby2_vx;
ILnew = IL;
IRnew = IR;

% resize the optical flow result to half
[h1,w1] = size(vx);
h2 = round(h1/2);
w2 = round(w1/2);
vx2 = zeros(h2,w2);
vy2 = vx2;

for i = 1:h2
    for j = 1:w2
      if 2*i<=h1 && 2*j<=w2
        t = vx(2*i-1:2*i,2*j-1:2*j);
        vx2(i,j) = mean(t(:));  
        t = vy(2*i-1:2*i,2*j-1:2*j);
        vy2(i,j) = mean(t(:));
      else
          
        
      end
     
    end
end



for i = 1:IH
    for j = 1:IW
     dx = round(vx(i,j));
     dy = round(vy(i,j));
     m = i+dy;
     n = j+dx;
     %   comparision by block
     r = 1;       
     if m-r>0 && n-r>0 && m+r<IH+1 && n+r<IW+1
%         s = I1(i,j);
%         t = I2(m,n);
%          if s>0.1&&t>0.1
%        if s*t>0
%             I1new(i,j) = max(s,t);
%             I2new(m,n) = max(s,t);
%        else 
%            I1new(i,j) = 0;     % if I1 is result of soft edge
%           % I2new(i,j) = 0;  
%        end
           block1 = IL(i-r:i+r,j-r:j+r);
           block2 = IR(m-r:m+r,n-r:n+r);
           t1 = mean(block1(:));
           t2 = mean(block2(:));
           ILnew(i,j) = max(t1,t2);
           ILnew(m,n) = max(t1,t2);



     end
    end
end

figure,imshow(ILnew);
figure,imshow(IRnew);
% imwrite(I1new,'bottle23_set0.png','png');


% for i = 1:IH
%     for j = 1:IW
%      dx = round(vx(i,j));
%      dy = round(vy(i,j));
%      m = i+dy;
%      n = j+dx;
%      if m>0&&n>0&&m<(IH+1)&&n<(IW+1)
%         s = I1new(i,j);
%         t = I2new(m,n);
%         I1new(i,j) = max(s,t);
%         I2new(m,n) = max(s,t);     
%      end
%        
%    end
% end
% figure,imshow(I1new);
% figure,imshow(I2new);
